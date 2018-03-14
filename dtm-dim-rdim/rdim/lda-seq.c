// Authors: David Blei (blei@cs.princeton.edu)
//          Sean Gerrish (sgerrish@cs.princeton.edu)
//          Yuening Hu (ynhu.moon@gmail.com)
//          Aaron Gerow (gerow@uchicago.edu)
//
// Copyright 2011 Sean Gerrish and David Blei
// All Rights Reserved.
//
// rDIM version: Copyright 2015 Aaron Gerow, Yuening Hu, Jordan Boyd-Graber, James Evans and David Blei
// All Rights Reserved.
//
// See the README for this package for details about modifying or
// distributing this software.

#include <unistd.h>
#include <errno.h>
#include "lda-seq.h"
#include "gflags.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "main.h"
#include "checkpoint.h"

#ifdef __cplusplus 
extern "C" { 
#endif 

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector_double.h>

#ifdef __cplusplus
} 
#endif

#include <omp.h>
#include <sys/types.h>
#include <sys/wait.h>

#ifdef __cplusplus 
extern "C" { 
#endif 
namespace mklspace {

#include "mkl_lapacke.h"
#include "mkl.h"
#include "mkl_types.h"
#include "mkl_spblas.h"
#include "mkl_pardiso.h"
}
#ifdef __cplusplus 
} 
#endif

#ifdef CUDA
#include <cuda_runtime.h>
#include <cublas.h>
#include <cusparse.h>
#endif

#ifdef SPARSE
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>
#include <gsl/gsl_spblas.h>
#endif

#define PREC 15
#define DGEMM_MKL 1
#define DGEMM_GSL 2
#define DGEMM_SPMKL 3
#define PREC 15

extern int     LDA_INFERENCE_MAX_ITER;
static int     DGEMM_METHOD=DGEMM_MKL;
static double* scaled_influence = NULL;
const  int     TIME = -3;
const  double  PI   =  3.14159265358979;

DEFINE_int32(lda_sequence_max_iter, 20,    "The maximum number of iterations.");
DEFINE_int32(lda_sequence_min_iter,  1,    "The minimum number of iterations.");
DEFINE_int32(forward_window,         1,    "The forward window for deltas." 
				           "If negative, we use a beta with mean 5.");
DEFINE_string(normalize_docs, "normalize", "Describes how documents's wordcounts " 
			                   "are considered for finding influence."
	      				   "Options are \"normalize\", \"none\", "
	      				   "\"occurrence\", \"log\", or \"log_norm\".");
DEFINE_int32(save_time,             -1,    "Save a specific time.  If -1, save all times.");
DEFINE_int32(fix_topics,             0,    "Fix a set of this many topics. This amounts " 
	                                   "to fixing these topics' variance at 1e-10.");
DEFINE_double(sigma_mu,             1.0,   "Hyper-prior defining the distribution, N(0,sigma_mu), from which to draw \\hat{\\mu}.");
DEFINE_double(conditioning_epsilon, 0.0,   "Amount of noise to add when inverting degenerate sparse matrices");
DEFINE_string(read_hdf_term1,       "",    "Path to the location of an HDF5 version of the ew32_term1");
DEFINE_string(read_csr_term1,       "",    "Path to the location of an CSR3 version of the ew32_term1");
DEFINE_bool(enable_checkpointing,   false, "If true, export an HDF5 checkpoint after every iteration.");
DEFINE_bool(die_after_checkpoint,   false, "Exit after the first checkpoint is written.");
DEFINE_string(checkpoint_outdir,    "",    "Write checkpoints to this directory.");
DEFINE_string(checkpoint_recover,   "",    "Start this run from the checkpoint stored here.");

DECLARE_string(model);
DECLARE_int64(rng_seed);
DECLARE_int32(max_number_time_points);
DECLARE_double(sigma_d);
DECLARE_double(sigma_l);
DECLARE_double(sigma_c);
DECLARE_double(sigma_cv);
DECLARE_int32(debug);
DECLARE_string(outname);
DECLARE_string(corpus_prefix);
DECLARE_int32(threads);
DECLARE_int32(kthreads);

#ifdef CUDA
void cusparse_errchk(const int err, const char* mssg) {
  if (err != CUSPARSE_STATUS_SUCCESS) {
    printf("!! cuSPARSE error (%s): ",mssg);

    switch (err) {
    case CUSPARSE_STATUS_NOT_INITIALIZED:
      printf("CUSPARSE_STATUS_NOT_INITIALIZED\n");
      break;
    case CUSPARSE_STATUS_ALLOC_FAILED:
      printf("CUSPARSE_STATUS_ALLOC_FAILED\n");
      break;
    case CUSPARSE_STATUS_INVALID_VALUE:
      printf("CUSPARSE_STATUS_INVALID_VALUE\n");
      break;
    case CUSPARSE_STATUS_ARCH_MISMATCH:
      printf("CUSPARSE_STATUS_ARCH_MISMATCH\n");
      break;
    case CUSPARSE_STATUS_MAPPING_ERROR:
      printf("CUSPARSE_STATUS_MAPPING_ERROR\n");
      break;
    case CUSPARSE_STATUS_EXECUTION_FAILED:
      printf("CUSPARSE_STATUS_EXECUTION_FAILED\n");
      break;
    case CUSPARSE_STATUS_INTERNAL_ERROR:
      printf("CUSPARSE_STATUS_INTERNAL_ERROR\n");
      break;
    case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
      printf("CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED\n");
    default :
      printf("Unknown Error...\n");
    }
  }
  
}
#endif

/*
 * populate an LDA model at a particular time point
 *
 */
inf_var* inf_var_alloc(int number_topics, corpus_seq_t* corpus_seq) {
  // Hate to do this, but I had trouble using it before.  This should
  // be the first place we use it; otherwise we'll get a sigsev nil. // SG

  if (scaled_influence == NULL)
    scaled_influence = NewScaledInfluence(FLAGS_max_number_time_points);

  inf_var* inf_var_ptr                  = (inf_var*)     malloc(sizeof(inf_var));
  inf_var_ptr->doc_weights              = (gsl_matrix**) malloc(sizeof(gsl_matrix*) * corpus_seq->len);
  inf_var_ptr->renormalized_doc_weights = (gsl_matrix**) malloc(sizeof(gsl_matrix*) * corpus_seq->len);
  inf_var_ptr->ntime                    = corpus_seq->len;
  int i=0;
  
  for (i=0; i < corpus_seq->len; ++i) {
    corpus_t* corpus = corpus_seq->corpus[i];
    outlog("creating matrix. %d %d\n", corpus->ndocs, number_topics);
    
    if (corpus->ndocs == 0) {
      inf_var_ptr->doc_weights[i]        = (gsl_matrix*) malloc(sizeof(gsl_matrix));
      inf_var_ptr->doc_weights[i]->size1 = 0;
      inf_var_ptr->doc_weights[i]->size2 = number_topics;
      inf_var_ptr->renormalized_doc_weights[i] = (gsl_matrix*) malloc(sizeof(gsl_matrix));
      inf_var_ptr->renormalized_doc_weights[i]->size1 = 0;
      inf_var_ptr->renormalized_doc_weights[i]->size2 = number_topics;
    } else {
      inf_var_ptr->doc_weights[i]              = gsl_matrix_calloc(corpus->ndocs, number_topics);
      inf_var_ptr->renormalized_doc_weights[i] = gsl_matrix_calloc(corpus->ndocs, number_topics);
    }
  }
  
  return inf_var_ptr;
}

// added by YH
void inf_var_free(inf_var* ptr) {
  for (int i=0; i < ptr->ntime; i++) {
    gsl_matrix_free(ptr->doc_weights[i]);
    gsl_matrix_free(ptr->renormalized_doc_weights[i]);
  }
  
  free(ptr->doc_weights);
  free(ptr->renormalized_doc_weights);
  free(ptr);
}

#ifdef SPARSE
#ifdef MKL

int MKLSolve_lapacke_single(gsl_matrix* gsl_A, const gsl_vector* gsl_b, gsl_vector* gsl_x) {
  MKL_INT
    n    = gsl_A->size1,
    nrhs = 1, info, ipiv[n];
  double rcond, berr[nrhs], rpivot, rpvgrw;
  char   equed = 'N';
  
  //  double err_bnds_norm[nrhs], ferr[nrhs], err_bnds_comp[nrhs];
  int n_err_bnds = 0, nparams = 1;
  double params[1] = {0.0}; // Disables iterative refinement because we only use this for small solves
  
  // Use this call to /enable/ iterative refinement (eg. ifndfef SPARSE calls this function
  // for a large solve on eq32_term1. Don't forget to set params[1] = 1.0 above
  /*  info = mklspace::LAPACKE_dgesvx(LAPACK_ROW_MAJOR, 'N', 'N',    n,  nrhs, gsl_A->data,  
				      n,     NULL,  n, ipiv, &equed, NULL, NULL,  gsl_b->data,  
				      nrhs , gsl_x->data, nrhs, &rcond, ferr, berr,   &rpivot);        
  */

  info = mklspace::LAPACKE_dgesvxx(LAPACK_ROW_MAJOR, 'N', 'N',    n,  nrhs, gsl_A->data,  
				   n,    NULL,  n, ipiv, &equed, NULL, NULL,  gsl_b->data,  
				   nrhs, gsl_x->data, nrhs, &rcond, &rpvgrw, berr, n_err_bnds,
				   NULL /*err_bnds_norm*/, NULL /*err_bnds_cmp*/, nparams, params);
 
  if (info == (n+1)) {
    fprintf(stderr, "(LAPACKE linear solve) WARNING: Condition number below MACHINE_EPSILON.\n", (int) info);
  }
  else if (info > 0) {
    fprintf(stderr, "(LAPACKE linear solve) Matrix looks degenerate (retcode %i).\n", (int) info);
  }
  else if (info < 0) {
    fprintf(stderr, "(LAPACKE linear solve) threw an error: %i, exiting...\n", (int) info);
    _exit(info);
  }

  return(info);
}

int MKLSolve_lapacke_comp(gsl_matrix* gsl_A, gsl_vector* gsl_b, gsl_vector* gsl_x) {
  MKL_INT 
    n    = gsl_A->size1, 
    nrhs = 1, 
    info, 
    ipiv[n];
  
  gsl_vector_memcpy(gsl_x, gsl_b);
  
  // Factor:
  info = mklspace::LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, gsl_A->data, n, ipiv);

  // Solve (result is set in gsl_x->data
  info = mklspace::LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', n, nrhs, gsl_A->data, n, ipiv, gsl_x->data, nrhs);

  if (info == (n+1))
    fprintf(stderr, "(LAPACKE linear solve) WARNING: Condition number below MACHINE_EPSILON.\n", (int) info);
  else if (info > 0)
    fprintf(stderr, "(LAPACKE linear solve) Matrix looks degenerate (retcode %i).\n", (int) info);
  else if (info < 0) {
    fprintf(stderr, "(LAPACKE linear solve) threw an error: %i, exiting...\n", (int) info);
    _exit(info);
  }

  return(info);
}

int MKLSolve_lapacke(gsl_matrix* gsl_A, gsl_matrix* gsl_b, gsl_matrix* gsl_x, 
		     mkl_factor* factor, bool& try_factor, bool& A_factored) {
  MKL_INT
    n    = gsl_A->size1,
    nrhs = gsl_b->size2, // b and x have to be COL_MAJOR -- despite arg[0] to dgesvx below
    info;
  double
    rcond, berr[nrhs], ferr[nrhs], rpivot, 
    epsilon = FLAGS_conditioning_epsilon ? 1.0 / pow(10, FLAGS_conditioning_epsilon) : 0.0;
  int    tries = 9;
  char   equed = 'N';
  time_t start;
  
  // A has not been factored (this is the first solve):
  if (!A_factored) {
    outlog("INFO: This appears to be the first time solving for mu: will try to save a factored form of eq32_term1, equilibrating scales and pivot tables for future use. However, this doubles the memory footprint, so if we fail here, that's probably why...%s\n","");
    
    if ((factor->af = gsl_matrix_alloc(n,n)) != NULL) { // The hard one
      factor->ipiv  = (MKL_INT*) malloc(sizeof(MKL_INT) * n);
      factor->r     = (double*)  malloc(sizeof(double)  * n);
      factor->c     = (double*)  malloc(sizeof(double)  * n);
      try_factor    = true;
    }
    else 
      outlog("BUMMER: Not enough space to cache the factorization from this solve, subsequent solves will take just as long. If you can ~double the memory for this run, it will go almost twice as fast.%s\n","");
  }
  else {
    outlog("INFO: Looks like we have a cached factorizaton from an earlier solve -- using that.%s\n","");
    try_factor = false;
  }
  
  do {
    outlog("INFO: Solving...\n%s","");
    
    if (A_factored) {
      info = mklspace::LAPACKE_dgesvx(LAPACK_ROW_MAJOR, 'F', 'N', n, nrhs, gsl_A->data, n, 
				      factor->af->data, n, factor->ipiv, &equed, factor->r, 
				      factor->c, gsl_b->data, nrhs, gsl_x->data, nrhs, 
				      &rcond, ferr, berr, &rpivot);
    }
    else if (try_factor) {
      debug("Timinng factor(A)::solve...%s\n","");
      info = mklspace::LAPACKE_dgesvx(LAPACK_ROW_MAJOR, 'N', 'N', n, nrhs, gsl_A->data, n, 
				      factor->af->data, n, factor->ipiv, &equed, factor->r, 
				      factor->c, gsl_b->data, nrhs, gsl_x->data, nrhs, 
				      &rcond, ferr, berr, &rpivot);
    }
    else {
      debug("Timinng non factor::solve...%s\n","");
      
      MKL_INT ipiv[n]; // Should be find on the stack ... famous last words
      info = mklspace::LAPACKE_dgesvx(LAPACK_ROW_MAJOR, 'N', 'N', n, nrhs, gsl_A->data, n, 
				      NULL, n, ipiv, &equed, NULL,
				      NULL, gsl_b->data, nrhs, gsl_x->data, nrhs, 
				      &rcond, ferr, berr, &rpivot);
    }
    
    outlog("...done (retcode = %i)\n", info);
    
    if (info == (n+1)) {
      outlog("LAPACKE WARDNING: Condition number below MACHINE_EPSILON (retcode %i).\n", info);
    }
    else if (info > 0) {
      outlog("LAPACKE WARNING: Matrix looks degenerate (retcode %i), adding noise...\n", info);
      add_noise_scaling(gsl_A, epsilon);
    }
    else if (info < 0) {
      outlog("(LAPACKE ERROR: %i (-n implies arg n to lapacke_dgesvx() is \"wrong\"), exiting...\n", info);
      _exit(info);
    }
  } while (--tries > 0 && info > 0 && info != (n+1));

  if (try_factor && info == 0) {
    A_factored = true;
    factor->n  = n;
  }

  debug("LAPACKE x[0-4]=%f %f %f %f\n",
	gsl_x->data[0],
	gsl_x->data[1],
	gsl_x->data[2],
	gsl_x->data[3]);

  return(info);
}

void hdf_write_matrix(const gsl_matrix* M, const char* fname) {
  hid_t fid;
  hsize_t dims[2] = {M->size1, M->size2};
  herr_t status;

  fid    = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  status = H5LTmake_dataset(fid, "term1", 2, dims, H5T_NATIVE_DOUBLE, M->data);
  H5LTset_attribute_ulong(fid, "term1", "size1", &M->size1, 1);
  H5LTset_attribute_ulong(fid, "term1", "size2", &M->size2, 1);
  assert(status != -1);

  status = H5Fclose(fid); assert(status != -1);
  assert(status != -1);
  return;
}

gsl_matrix* read_hdf_term1(const int s) {
  outlog("Reading HDF5 term 1...%s\n","")
  
  char fname[BUFSIZ];
  sprintf(fname, "%s", FLAGS_read_hdf_term1.c_str());
  
  gsl_matrix* ret = gsl_matrix_alloc(s, s);
  hid_t       fid;
  herr_t      status;

  fid    = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
  status = H5LTread_dataset_double(fid, "term1", ret->data);
  
  status = H5Fclose(fid);
  if (status == -1) {
    outlog("HDF5 returned -1 on reading: %i\n", (int) status);
    _exit(1);
  }

  int NZ = 0;
  for (long i=0; i<s*s; i++)
    if (ret->data[i] != 0.0)
      NZ++;

  debug("read in (%i x %i) matrix with NZ=%i\n",ret->size1,ret->size2, NZ);
  return ret;
}

gsl_matrix* hdf_read_matrix(const char* fname) {
  gsl_matrix*   ret;
  hid_t         fid;
  herr_t        status;
  unsigned long s1, s2;
  
  fid    = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
  H5LTget_attribute_ulong(fid, "term1", "size1", &s1);
  H5LTget_attribute_ulong(fid, "term1", "size2", &s2);

  if (s1 < 1 || s2 < 1) {
    fprintf(stderr,"ERROR: Got matrix->size1 < 1 for %s (%i x %i), exiting...\n", fname, s1,s2);
    _exit(24);
  }

  ret = gsl_matrix_alloc(s1, s2);

  status = H5LTread_dataset_double(fid, "term1", ret->data);
  status = H5Fclose(fid);

  if (status == -1) {
    fprintf(stderr,"HDF5 returned -1 on reading %s: %i\n", fname, (int) status);
    _exit(22);
  }

  return ret;
}

gsl_matrix* read_hdf_b(const int bs1, const int bs2, const char* fname) {
  outlog("Reading HDF5 term b...%s\n","");
  
  gsl_matrix* ret = gsl_matrix_alloc(bs1,bs2);
  hid_t       fid;
  herr_t      status;

  fid    = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
  H5LTread_dataset_double(fid, "term1", ret->data);
  
  status = H5Fclose(fid);
  if (status == -1) {
    outlog("HDF5 returned -1 on reading: %i\n", (int) status);
    _exit(1);
  }

  return ret;
}

gsl_matrix* read_hdf_A(const int s, const char* fname) {
  outlog("Reading HDF5 term A...%s\n","");
  
  gsl_matrix* ret = gsl_matrix_alloc(s,s);
  hid_t       fid;
  herr_t      status;

  fid    = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
  H5LTread_dataset_double(fid, "term1", ret->data);
  
  status = H5Fclose(fid);
  if (status == -1) {
    outlog("HDF5 returned -1 on reading: %i\n", (int) status);
    _exit(1);
  }

  return ret;
}

int hdf_write_csr3(mkl_spmat* S, const char* fname) {
  hid_t fid;
  hsize_t  dims[1] = {S->NZ};
  hsize_t idims[1] = {(S->n)+1};
  herr_t status;
  MKL_INT n1 = (S->n)+1;

  debug("WRITE S2: In hdf_write_csr3\nNZ=%i,ja[1919]=%i ja[0]=%i n1=%i\n",
	S->NZ, S->Mj[1919], S->Mj[0], n1);

  for (long i = 0; i < S->NZ; i++) {
    if (S->Mj[i] < 0)
      debug("ERROR SS1: Sj[%i]=%i (v=%f)\n",i,S->Mj[i],S->Mv[i]);
  }

  fid    = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  status = H5LTmake_dataset(  fid, "CSR3_a",  1,   dims, H5T_NATIVE_DOUBLE, S->Mv);
  status = H5LTmake_dataset(  fid, "CSR3_ia", 1,  idims, H5T_NATIVE_LLONG,  S->Mi);
  status = H5LTmake_dataset(  fid, "CSR3_ja", 1,   dims, H5T_NATIVE_LLONG,  S->Mj);
  H5LTset_attribute_long_long(fid, "CSR3_a", "NZ",  &S->NZ, 1);
  H5LTset_attribute_long_long(fid, "CSR3_a", "n1",  &n1, 1);
  assert(status != -1);
  
  status = H5Fclose(fid);
  assert(status != -1);
  
  return (int) status;

}

int hdf_write_gsl_to_csr3(gsl_matrix* M, const char* fname) {
  MKL_INT 
    n  = (MKL_INT) M->size1,
    m  = (MKL_INT) M->size2,
    NZ = 0,     
    c  = 0, 
    jc = 0;
  double v = 0;

    for (MKL_INT i=0; i<n; i++)
      for (MKL_INT j=i; j<m; j++) // Only upper triangle
	if (gsl_matrix_get(M,i,j) != 0.0) 
	  NZ++;
  debug("In CSR3 writer...NZ=%i, fname=%s\n",NZ,fname);

  MKL_INT
    *ia = (MKL_INT*) malloc(sizeof(MKL_INT) * (n+1)),
    *ja = (MKL_INT*) malloc(sizeof(MKL_INT) * NZ);
  double 
    *a  = (double*) malloc(sizeof(double) * NZ); 

  debug("Writing...1...%s","");
  bool rowset;
  for (MKL_INT i=0; i<n; i++) {
    rowset = false;
    for (MKL_INT j=i; j<m; j++) { // Note this is only the upper triangle
      v = gsl_matrix_get(M,i,j);

      if (v!=0.0) {
	a[c]=v;
	ja[c]=j;

	if (!rowset) {
	  ia[jc++] = c;
	  rowset=true;
	}

	c++;
      }
    }
  }
  ia[n]=NZ;
  debug("2...%s","");

  debug("a%s","");
  hsize_t  a_dims[1] = {NZ};
  debug("b%s","");
  hsize_t ia_dims[1] = {n+1};
  debug("c%s","");
  herr_t status;
  debug("d%s","");

  MKL_INT n1 = n+1;
  debug("e%s","");
  debug("f%s","");
  hid_t fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  debug("3...%s","");
  debug("b%s","");
  status = H5LTmake_dataset(fid, "CSR3_a",  1,  a_dims, H5T_NATIVE_DOUBLE, a);
  debug("4...%s","");
  status = H5LTmake_dataset(fid, "CSR3_ia", 1, ia_dims, H5T_NATIVE_LLONG, ia);
  debug("5...%s","");
  status = H5LTmake_dataset(fid, "CSR3_ja", 1,  a_dims, H5T_NATIVE_LLONG, ja);
  debug("6...%s","");
  H5LTset_attribute_long_long(  fid, "CSR3_a", "NZ", &NZ, 1);
  debug("a%s","");
  H5LTset_attribute_long_long(  fid, "CSR3_a", "n1", &n1, 1);
  debug("b%s","");
  assert(status != -1);
  debug("7...%s","");  
  status = H5Fclose(fid);
  debug("a%s","");
  assert(status != -1);
  debug("done!\n%s","");
  return (int) status;
}

int fork_pardiso_solve(const char* A_fname, gsl_matrix* &b, gsl_matrix* &x, const int cNZ) {
  char fname[BUFSIZ];
  char wd[BUFSIZ];
  char Aloc[BUFSIZ];
  sprintf(wd, "--wd=%s",FLAGS_outname.c_str());
  sprintf(Aloc, "--Aloc=%s",A_fname);

  sprintf(fname, "%s/b.h5",FLAGS_outname.c_str());  
  hdf_write_matrix(b,fname);

  pid_t parent = getpid();
  outlog("Forking to MKL::PARDIISO Solver%s\n\n","");

  pid_t pid = fork();

  if (pid == -1) {
    outlog("FAILURE: fork() to MKL::PARDISO failed...%s\n","");
    return -1;
  } 
  else if (pid > 0)  {
    int status;
    waitpid(pid, &status, 0);
  }
  else {
    // Have to do some seriouly MAD threading magic here so the solver has access tothreads 
    // exposed by libgomp. If compiled with icc or clang, this mightn't work...
    cpu_set_t *mask;
    size_t size;
    int nrcpus = 256; 

    mask = CPU_ALLOC(nrcpus);
    size = CPU_ALLOC_SIZE(nrcpus);
    for (int i = 0; i < nrcpus; i++)
      CPU_SET_S(i, size, mask);

    if (sched_setaffinity(0, size, mask) == -1) { 
      outlog("ERROR: This is forked main (pre-exec()): Failed to to preserve GOMP thread-mask. Big problem..., can't even exit here.%s\n","");
      _exit(EXIT_FAILURE);
    }
    
    CPU_FREE(mask);
    execlp("./pardiso_solver", "./pardiso_solver", wd, Aloc, (char*)0);
    _exit(EXIT_FAILURE); // exec should not return anything here.
  }  
  
  outlog("*** This is main (DTM/DIM/rDIM model) ***%s\n","");
  sprintf(fname, "%s/x.h5",FLAGS_outname.c_str());  
  x = hdf_read_matrix(fname);
  
  return 0;
}

void write_gsl_matrix_col_major(const gsl_matrix* M, const char* fname) {
  FILE* f = fopen(fname, "w");
  
  for (int j = 0; j < M->size2; j++ )
    for (int i = 0; i < M->size1; i++ ) 
      fprintf(f, "%.14f\n", gsl_matrix_get(M,i,j));

  fclose(f);
}

#endif

int SparseSolve(const gsl_spmatrix* A, const gsl_vector* b, gsl_vector* x, const int s) {
  const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
  size_t iter;
  bool report = true;
  int
    m,
    max_super_tries = 5,
    max_iter  = 25, // Number of time to iterate the GMRES solver
    status    = GSL_CONTINUE;
  double 
    residual, 
    tol     = 1.0e-16;
  
  while (max_super_tries-- >= 0 && status != GSL_SUCCESS) {
    m             = s;
    int max_tries = 10; // Each try will increase max_iter and m
    max_iter      = 25, // Number of time to iterate the GMRES solver
    
    assert((int) pow(max_tries, m)); // This should not exceed limits.h::INT_MAX on a given system.
    
    do {
      iter = 0;
      gsl_splinalg_itersolve *work = gsl_splinalg_itersolve_alloc(T, b->size, m);
      
      do
        status = gsl_splinalg_itersolve_iterate(A, b, tol, x, work);
      while (status == GSL_CONTINUE && ++iter < max_iter);
  
      gsl_splinalg_itersolve_free(work);
  
      if (status != GSL_SUCCESS) {
        outlog("WARNING: SparseSolve() failed to converge in %d iterations (sub-space m=%d)...increasing iterations and sub-space dimensions.\n", max_iter, m);
        m        *= 1.5;
	m         = m >= b->size ? (int) b->size : m;
        max_iter *= 1.5;
        report    = true; // Because we failed once, users might want to know when (if) we succeed.
      }
    } while (status != GSL_SUCCESS && max_tries > 0 && max_tries--);
    if (status != GSL_SUCCESS) {
      outlog("RECOVERABLE FAILURE: SparseSolve() failed to converge in %d iterations (sub-space m=%d), lowering tolerance to %g.\n", max_iter, m, tol);
      tol *= 10;
    }
    else if (report) {
      outlog("SUCCESS: SparseSolve() converged in %d iterations (sub-space m=%d).\n", (int) iter, m);
    }
  }

  return status;
}
#endif

void LUSolve(gsl_matrix* A, const gsl_vector* b, gsl_vector* x) {
  int permutation_sign;
  gsl_permutation* permutation = gsl_permutation_alloc(b->size);
  
  gsl_linalg_LU_decomp(A, permutation, &permutation_sign);
  gsl_linalg_LU_solve(A, permutation, b, x);
  
  gsl_permutation_free(permutation);
}

// Solves the linear system Ax = b for x.
// Assumes that x is already allocated.
void CholeskySolve(gsl_matrix* A, const gsl_vector* b, gsl_vector* x) {
  gsl_linalg_cholesky_decomp(A);
  gsl_linalg_cholesky_solve(A, b, x);
}

// Find the sums of influence of all documents in advance.
// g has scores for everything *up to but not including* g.
// AG: What is this doing and why?
void InfluenceTotalFixed(lda_seq* seq) {
  gsl_vector* exp_tmp = gsl_vector_alloc(seq->nterms);

  for (int k = 0; k < seq->ntopics; ++k) { // For every topic
    for (int s = 0; s < seq->nseq; ++s) {  // For every time-slice
      gsl_vector_view g = gsl_matrix_column(seq->influence_sum_lgl[k], s);
      gsl_vector_set_zero(&g.vector);
      
      for (int t=0; t <= s; ++t) { // For every slice, t, up to current slice s.
	gsl_vector_view w_phi_l = gsl_matrix_column(seq->topic[k]->w_phi_l, t);
	gsl_vector_memcpy(exp_tmp, &w_phi_l.vector);
	gsl_vector_scale( exp_tmp,  scaled_influence[s - t]);
	gsl_vector_add(  &g.vector, exp_tmp);
      }
    }
  }
  
  gsl_vector_free(exp_tmp);
}

void DumpTimeDocTopicStats(const char* root, size_t t, corpus_t* corpus, gsl_matrix** phi) {
  char name[400];
  sprintf(name, "%s%ld_doc_term_topics.dat", root, t);
  FILE* f = fopen(name, "w");

  for (unsigned int d=0; d < corpus->ndocs; ++d) {
    gsl_matrix* phi_d = phi[d];
    doc_t* doc        = corpus->doc[d];

    for (unsigned int n=0; n < doc->nterms; ++n) {
      unsigned int w = doc->word[n];

      // First, find the max topic weight.
      unsigned int max_topic_index  = 0;
      double       max_topic_weight = gsl_matrix_get(phi_d, n, 0);

      for (unsigned int k=0; k < phi_d->size2; ++k) {
	double phi_d_n_k = gsl_matrix_get(phi_d, n, k);

	if (phi_d_n_k > max_topic_weight) {
	  max_topic_weight = phi_d_n_k;
	  max_topic_index  = k;
	}
      }
      
      fprintf(f, "%d:%d:%.3f", w, max_topic_index, max_topic_weight);
      
      if (n < doc->nterms - 1)
	fprintf(f, " ");
    }

    fprintf(f, "\n");
  }

  fclose(f);
}

#ifdef CUDA
int cublas_dgemm(gsl_matrix *A, const gsl_matrix* B, gsl_matrix* C) {
  // Target is C=A^T*B
  assert(A->size2 == C->size1);
  assert(B->size2 == C->size1);
  assert(A->size1 == B->size1);
  
  const size_t NA     = A->size1;
  const size_t MA     = A->size2;
  const size_t NB     = B->size1;
  const size_t MB     = B->size2;
  const size_t NC     = C->size1;
  const size_t MC     = C->size2;
        int    status = NULL;
	double alpha  = 1.0;
	double beta   = 0.0;
	double *cA, *cB, *cC;
	double *CT    = (double*) malloc(sizeof(double) * (NC*MC));

	/*
	  status = cublasInit();
  if (status != CUBLAS_STATUS_SUCCESS) {
      fprintf(stderr, "!!!! CUBLAS initialization error\n");
      
      return EXIT_FAILURE;
  }
	*/
  status = cudaMalloc(&cA, NA*MA * sizeof(*cA));
  //  status = cublasAlloc(NA*MA, sizeof(*cA), (void**) &cA);
  if (status != CUBLAS_STATUS_SUCCESS) {
      fprintf (stderr, "!!!! device memory allocation error (A)\n");
      return EXIT_FAILURE;
  }
 
  status = cudaMalloc(&cB, NB*MB * sizeof(*cB));
  if (status != CUBLAS_STATUS_SUCCESS) {
      fprintf (stderr, "!!!! device memory allocation error (B)\n");
      return EXIT_FAILURE;
  }
 
  status = cudaMalloc(&cC, NC*MC * sizeof(*cC));
  if (status != CUBLAS_STATUS_SUCCESS) {
      fprintf (stderr, "!!!! device memory allocation error (C)\n");
      return EXIT_FAILURE;
  }
  
  /* Initialize the device matrices with the host matrices */
  status = cublasSetVector(NA*MA, sizeof(*cA), A->data, 1, cA, 1);
  if (status != CUBLAS_STATUS_SUCCESS) {
      fprintf (stderr, "!!!! device access error (write A)\n");
      return EXIT_FAILURE;
  }
 
  status = cublasSetVector(NB*MB, sizeof(*cB), B->data, 1, cB, 1);
  if (status != CUBLAS_STATUS_SUCCESS) {
      fprintf (stderr, "!!!! device access error (write B)\n");
      return EXIT_FAILURE;
  }

  cublasGetError();
   
  // cuBLAS is col-major, so transpose B on multiply, later we transpose C...
  cublasDgemm('n', 't',
	      NC, MC, NB,
	      alpha, cA, MA,
	             cB, MB, beta,
	             cC, NC);
  
  status = cublasGetError();
  
  if (status != CUBLAS_STATUS_SUCCESS) {
    plainlog("!!!! kernel execution error.\n");
    outlog("!!!! errno %d.\n", status);
    return EXIT_FAILURE;
  }
 
  status = cublasGetVector(NC*MC, sizeof(*CT), cC, 1, CT, 1);
  if (status != CUBLAS_STATUS_SUCCESS) {
      fprintf (stderr, "!!!! device access error (read C)\n");
      return EXIT_FAILURE;
  }
  
  status = cudaFree(cA);
  if (status != CUBLAS_STATUS_SUCCESS) {
      fprintf (stderr, "!!!! device access error (free A)\n");
      return EXIT_FAILURE;
  }

  status = cudaFree(cB);
  if (status != CUBLAS_STATUS_SUCCESS) {
      fprintf (stderr, "!!!! device access error (free B)\n");
      return EXIT_FAILURE;
  }
  
  status = cudaFree(cC);
  if (status != CUBLAS_STATUS_SUCCESS) {
      fprintf (stderr, "!!!! device access error (free C)\n");
      return EXIT_FAILURE;
  }

  double* c = CT;
  for (int i = 0 ; i < C->size1; i++)
    for (int j = 0; j < C->size2; j++)
      C->data[i+(j*C->size1)] = *c++; //CT[(i*N)+j];
    
  free(CT);

  //  cublasShutdown();
  
  return 1;
}

int cusparse_dgemm1(gsl_matrix *A, const gsl_matrix* B, gsl_matrix* C) {
  assert(A->size2 == C->size1);
  assert(B->size2 == C->size1);
  assert(A->size1 == B->size1);
  
  const int  an    = A->size1;
  const int  am    = A->size2;
  const int  bn    = B->size1;
  const int  bm    = B->size2;
  const int  cn    = C->size1;
  const int  cm    = C->size2;
        int  ldc   = C->size1; 

  const double alpha = 1.0;
  const double beta  = 0.0;
  double *CT   = (double*) malloc(sizeof(double) * (cn*cm));

  cusparseOperation_t notrans = CUSPARSE_OPERATION_NON_TRANSPOSE;
  cusparseOperation_t trans   = CUSPARSE_OPERATION_TRANSPOSE;

  cusparseStatus_t stat;
  //  cusparseHandle_t hndl;
  cudaError_t      err;
  char            *errs;
  cusparseMatDescr_t descrA;
  
  int     *gpu_AI, *gpu_AJ, *gpu_rownnzA;
  double  *gpu_Acsr, *gpu_A,*gpu_B, *gpu_C;
  
  int nnzA = 0;
  int *rownnzA = (int*) malloc(sizeof(int) * am);

  for (size_t j = 0; j < am; j++) {
    rownnzA[j] = 0;
    for (size_t i = 0; i < an; i++)
      if (gsl_matrix_get(A,i,j))
	rownnzA[j]++;

    nnzA += rownnzA[j];
  }
  
  // set cusparse matrix types
  //  stat = cusparseCreate(&hndl);
  cusparse_errchk(stat,"init");

  cudaMalloc(&gpu_rownnzA,  am*sizeof(int));
  cudaMalloc(&gpu_A,     an*am*sizeof(double));
  cudaMalloc(&gpu_B,     bn*bm*sizeof(double));
  cudaMalloc(&gpu_C,     cn*cm*sizeof(double));
  cudaMalloc(&gpu_Acsr,   nnzA*sizeof(double));
  cudaMalloc(&gpu_AI,   (am+1)*sizeof(int));
  cudaMalloc(&gpu_AJ,     nnzA*sizeof(int));
  
  err = cudaGetLastError();
  if (err != cudaSuccess) { printf("!! [1] cudaError='%s'\n",cudaGetErrorString(err)); }
  
  cudaMemcpy(gpu_A,   A->data,   (an*am)*sizeof(double),    cudaMemcpyHostToDevice);
  err = cudaGetLastError();
  if (err != cudaSuccess) { printf("!! [2.0] cudaError='%s'\n",cudaGetErrorString(err)); }
  
  cudaMemcpy(gpu_rownnzA,   rownnzA,   (am)*sizeof(int),    cudaMemcpyHostToDevice);
  err = cudaGetLastError();
  if (err != cudaSuccess) { printf("!! [2.1] cudaError='%s'\n",cudaGetErrorString(err)); }

  cudaMemcpy(gpu_B,   B->data,   (bn*bm)*sizeof(double),    cudaMemcpyHostToDevice);
  err = cudaGetLastError();
  if (err != cudaSuccess) { printf("!! [2.2] cudaError='%s'\n",cudaGetErrorString(err)); }
  
  stat = cusparseCreateMatDescr(&descrA);
  cusparse_errchk(stat,"alloc");
  stat = cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparse_errchk(stat,"typing");
  stat = cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);
  cusparse_errchk(stat,"basing");

  err = cudaGetLastError();
  if (err != cudaSuccess) { printf("!! [3] cudaError='%s'\n",cudaGetErrorString(err)); }
  
  plainlog("7.2 Converting to CSR: A...");
    
  err = cudaGetLastError();
  if (err != cudaSuccess) { printf("!! [4] cudaError='%s'\n",cudaGetErrorString(err)); }

  stat = cusparseDdense2csr(hndl,am,an,descrA,gpu_A,am,gpu_rownnzA, // Inputs
			    gpu_Acsr, gpu_AI, gpu_AJ);              // Outputs

  err = cudaGetLastError();
  if (err != cudaSuccess) { printf("!! [5] cudaError='%s'\n",cudaGetErrorString(err)); }
  cusparse_errchk(stat,"dense2csr(A)");
      
  plainlog("\n7.3 Multiplying...\n");
  stat =  cusparseDcsrmm2(hndl, notrans, trans, cn, cm, bn, nnzA, &alpha,
			  descrA, gpu_Acsr, gpu_AI, gpu_AJ,
			  gpu_B, bm, &beta, gpu_C, cn);
  err = cudaGetLastError();
  if (err != cudaSuccess) { printf("!! [6] cudaError='%s'\n",cudaGetErrorString(err)); }
  cusparse_errchk(stat,"csrmm.");

  plainlog("7.4 Copying back to host...\n");

  cudaMemcpy(CT, gpu_C, cn*cm*sizeof(double), cudaMemcpyDeviceToHost);
  cusparse_errchk(stat,"Copy(A) from VRAM to main RAM.");

  plainlog("7.5 Transposing result...\n");

  double* c = CT;
  for (size_t i = 0 ; i < C->size1; i++)
    for (size_t j = 0; j < C->size2; j++)
      C->data[i+(j*C->size1)] = *c++; //CT[(i*N)+j];
  
  cudaFree(gpu_A);
  cudaFree(gpu_B);
  cudaFree(gpu_C);
  cudaFree(gpu_AI);
  cudaFree(gpu_AJ);
  cudaFree(gpu_Acsr);
  cudaFree(gpu_rownnzA);
  cusparseDestroyMatDescr(descrA);
  //  cusparseDestroy(hndl);
  
  free(CT);
  free(rownnzA);
  
  return 1;
}

// PUT cuSPARSE dgemm2 here!
int cusparse_dgemm2(gsl_matrix *A, const gsl_matrix* B, gsl_matrix* C) {
  assert(A->size2 == C->size1);
  assert(B->size2 == C->size1);
  assert(A->size1 == B->size1);
  
  const int  an    = A->size1;
  const int  am    = A->size2;
  const int  bn    = B->size1;
  const int  bm    = B->size2;
  const int  cn    = C->size1;
  const int  cm    = C->size2;
  double    *CT    = (double*) malloc(sizeof(double) * (cn*cm));

  cusparseOperation_t notrans = CUSPARSE_OPERATION_NON_TRANSPOSE;
  cusparseOperation_t trans   = CUSPARSE_OPERATION_TRANSPOSE;
  cusparseStatus_t    stat;
  //  cusparseHandle_t    hndl;
  cudaError_t         err;
  char               *errs;
  cusparseMatDescr_t  descrA, descrB, descrC;
 
  int     *gpu_AI, *gpu_BI,
          *gpu_AJ, *gpu_BJ,
          *gpu_CI, *gpu_CJ;
  double  *gpu_Acsr, *gpu_Bcsr, *gpu_Ccsr,
          *gpu_A,    *gpu_B, *gpu_C;
  
  int nnzA = 0, nnzB = 0, nnzC, baseC;
  int *gpu_rownnzA, *gpu_rownnzB, *gpu_rownnzC;
  int
    *rownnzA            = (int*) malloc(sizeof(int) * am),
    *rownnzB            = (int*) malloc(sizeof(int) * bm),
    *nnzTotalDevHostPtr = &nnzC;
  
  for (size_t j = 0; j < am; j++) {
    rownnzA[j] = 0;
    for (size_t i = 0; i < an; i++)
      if (gsl_matrix_get(A,i,j))
	rownnzA[j]++;

    nnzA += rownnzA[j];
  }
  
  // The transpoing count...  
  for (size_t j = 0; j < bm; j++) {
    rownnzB[j] = 0;
    for (size_t i = 0; i < bn; i++)
      if (gsl_matrix_get(B,i,j))
	rownnzB[j]++;

    nnzB += rownnzB[j];
  }

  //  stat = cusparseCreate(&hndl);

  cudaMalloc(&gpu_rownnzA,  am*sizeof(int));
  cudaMalloc(&gpu_rownnzB,  bm*sizeof(int));
  cudaMalloc(&gpu_A,     an*am*sizeof(double));
  cudaMalloc(&gpu_B,     bn*bm*sizeof(double));
  cudaMalloc(&gpu_C,     cn*cm*sizeof(double));
  cudaMalloc(&gpu_Acsr,   nnzA*sizeof(double));
  cudaMalloc(&gpu_Bcsr,   nnzB*sizeof(double));
  cudaMalloc(&gpu_AI,   (am+1)*sizeof(int));
  cudaMalloc(&gpu_BI,   (bm+1)*sizeof(int));
  cudaMalloc(&gpu_CI,   (cm+1)*sizeof(int));
  cudaMalloc(&gpu_AJ,     nnzA*sizeof(int));
  cudaMalloc(&gpu_BJ,     nnzB*sizeof(int));
    
  cudaMemcpy(gpu_A,       A->data,   (an*am)*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_rownnzA, rownnzA,      (am)*sizeof(int),    cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_B,       B->data,   (bn*bm)*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_rownnzB, rownnzB,      (bm)*sizeof(int),    cudaMemcpyHostToDevice);
  
  stat = cusparseCreateMatDescr(&descrA);
  stat = cusparseCreateMatDescr(&descrB);
  stat = cusparseCreateMatDescr(&descrC);

  stat = cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
  stat = cusparseSetMatType(descrB, CUSPARSE_MATRIX_TYPE_GENERAL);
  stat = cusparseSetMatType(descrC, CUSPARSE_MATRIX_TYPE_GENERAL);

  stat = cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);
  stat = cusparseSetMatIndexBase(descrB, CUSPARSE_INDEX_BASE_ZERO);
  stat = cusparseSetMatIndexBase(descrC, CUSPARSE_INDEX_BASE_ZERO);
  
  //  plainlog("8.2 Converting to CSR: A...");
    
  stat = cusparseDdense2csr(hndl,am,an,descrA,gpu_A,am,gpu_rownnzA, // Inputs
			    gpu_Acsr, gpu_AI, gpu_AJ);              // Outputs
  
  //  plainlog("B...\n");
  stat = cusparseDdense2csr(hndl,bm,bn,descrB,gpu_B,bm,gpu_rownnzB, // Inputs
			    gpu_Bcsr, gpu_BI, gpu_BJ);              // Outputs

  //  plainlog("8.3 Getting sparsity pattern of resultant...\n");
  
  cusparseSetPointerMode(hndl, CUSPARSE_POINTER_MODE_HOST);
  stat = cusparseXcsrgemmNnz(hndl,   notrans, trans, am, bm, cn,
			     descrA, nnzA, gpu_AI, gpu_AJ,
			     descrB, nnzB, gpu_BI, gpu_BJ,
			     descrC, gpu_CI, nnzTotalDevHostPtr);
  
  if (NULL != nnzTotalDevHostPtr)
    nnzC = *nnzTotalDevHostPtr;
  else {
    cudaMemcpy(&nnzC,  gpu_CI+cm, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&baseC, gpu_CI,    sizeof(int), cudaMemcpyDeviceToHost);
    nnzC -= baseC;
  }
  
  cudaMalloc((void**)&gpu_CJ,   sizeof(int)   *nnzC);
  cudaMalloc((void**)&gpu_Ccsr, sizeof(double)*nnzC);
  
  //  plainlog("8.4 Multiplying...\n");
  stat = cusparseDcsrgemm(hndl,   notrans, trans, am, bm, cn,
			  descrA, nnzA, gpu_Acsr, gpu_AI, gpu_AJ,
			  descrB, nnzB, gpu_Bcsr, gpu_BI, gpu_BJ,
			  descrC, gpu_Ccsr, gpu_CI, gpu_CJ);
  
  //  plainlog("8.5 Densifying result...\n");
  stat = cusparseDcsr2dense(hndl, cn, cm,
			    descrC, gpu_Ccsr, gpu_CI, gpu_CJ,
			    gpu_C, cn);

  //plainlog("8.6 Copying back to host...\n");

  cudaMemcpy(CT, gpu_C, cn*cm*sizeof(double), cudaMemcpyDeviceToHost);

  double* c = CT;
  for (size_t i = 0 ; i < C->size1; i++) {
    for (size_t j = 0; j < C->size2; j++) {
      gsl_matrix_set(C,j,i,*c++);
    }
  }
  
  cudaFree(gpu_rownnzA);
  cudaFree(gpu_rownnzB);
  cudaFree(gpu_A);
  cudaFree(gpu_B);
  cudaFree(gpu_C);
  cudaFree(gpu_Acsr);
  cudaFree(gpu_Bcsr);
  cudaFree(gpu_AI);
  cudaFree(gpu_BI);
  cudaFree(gpu_CI);
  cudaFree(gpu_AJ);
  cudaFree(gpu_BJ);
  cudaFree(gpu_CJ);
  cudaFree(gpu_Ccsr);
  
  cusparseDestroyMatDescr(descrA);
  cusparseDestroyMatDescr(descrB);
  cusparseDestroyMatDescr(descrC);
  //  cusparseDestroy(hndl);
  
  free(CT);
  free(rownnzA);
  free(rownnzB);
  
  return 1;
}

#endif

#ifdef MKL
int mkl_cblas_dgemm(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB,
                         double alpha, const gsl_matrix * A, const gsl_matrix * B,
                         double beta, gsl_matrix * C) {
  const size_t M      = C->size1;
  const size_t N      = C->size2;
  const size_t MA     = (TransA == CblasNoTrans) ? A->size1 : A->size2;
  const size_t NA     = (TransA == CblasNoTrans) ? A->size2 : A->size1;
  const size_t MB     = (TransB == CblasNoTrans) ? B->size1 : B->size2;
  const size_t NB     = (TransB == CblasNoTrans) ? B->size2 : B->size1;
  const int    transa = (TransA == CblasNoTrans) ? 111      : 112;
  const int    transb = (TransB == CblasNoTrans) ? 111      : 112;

  if (M == MA && N == NB && NA == MB) {
    mklspace::cblas_dgemm(mklspace::CblasRowMajor,
                          (mklspace::CBLAS_TRANSPOSE) transa,
                          (mklspace::CBLAS_TRANSPOSE) transb,
                          int (M), int (N), int (NA), alpha,
                          A->data, int (A->tda), B->data,
                          int (B->tda), beta, C->data, int (C->tda));
    return GSL_SUCCESS;
  }
  else {
    GSL_ERROR ("invalid length", GSL_EBADLEN);
  }
}

void mkl_spblas_dgemm(gsl_matrix *A, const gsl_matrix* B, gsl_matrix* C) {
  debug("SPMKL 1%s",""); 

  const MKL_INT an = A->size1;
  const MKL_INT am = A->size2;
  const MKL_INT bn = B->size1;
  const MKL_INT bm = B->size2;
  MKL_INT      ldc = C->size1; 
  MKL_INT       NZ = 0;
  MKL_INT*      AI = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * (an+1),64);
  MKL_INT*      BI = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * (bn+1),64);
  MKL_INT     info;       // Output status from MKL
  const char trans = 'T'; // Transpose the 'A' argument to MKL. 
                          // We actually transpose B because the col- to row-major 
                          // translation transposes A and C. Remember, this has to be C:=A^T * B

  // Note that this is not a general solver because NZ(A) is expected to equal NZ(B):
  double* Acsr;
  double* Bcsr;
  MKL_INT*  AJ;
  MKL_INT*  BJ;
  
  debug(", 2%s",""); 
#pragma omp parallel shared(NZ)
  {
#pragma omp for reduction(+:NZ) schedule(static)
    for (long i = 0; i < an*am; i++)
      if (A->data[i] != 0.0)
	NZ = NZ+1;
  }
  
  debug(", 3%s","");
  MKL_INT job[6] = {0, 0, 1, 2, NZ, 1};
  Acsr   = (double*)  mklspace::mkl_malloc(sizeof(double)  * NZ, 64);
  Bcsr   = (double*)  mklspace::mkl_malloc(sizeof(double)  * NZ, 64);
  AJ     = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * NZ, 64);
  BJ     = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * NZ, 64);
  
  debug(", 4%s","");

#pragma omp parallel num_threads(2)
  {
#pragma omp task
  mklspace::mkl_ddnscsr(job, &an, &am, A->data, &am, Acsr, AJ, AI, &info);
#pragma omp task
  mklspace::mkl_ddnscsr(job, &bn, &bm, B->data, &bm, Bcsr, BJ, BI, &info);
  }

  debug(", 5%s","");
  //omp_set_num_threads(1);
  //omp_set_dynamic(0);
  //omp_set_dynamic(1); // Old way had this uncommented
  debug("a%s","");
  //mklspace::mkl_set_dynamic(1); // Old way had this uncommented
  //mklspace::mkl_set_num_threads(FLAGS_threads);
  //mklspace::mkl_set_num_threads_local(FLAGS_threads);
  //omp_set_num_threads(FLAGS_threads);
  debug("b%s","");
  debug(" [GOMP=%i, MKL=%i] ", omp_get_num_threads(), mklspace::mkl_get_max_threads());

  mklspace::mkl_dcsrmultd(&trans, &bn, &bm, &am, Bcsr, BJ, BI, Acsr, AJ, AI, C->data, &ldc);
  debug("c%s","");
  //omp_set_dynamic(1);
  debug("d%s","");
  //mklspace::mkl_set_num_threads(1);
  //omp_set_num_threads(FLAGS_threads);

  debug(", 6 %s","");
  mklspace::mkl_free(Acsr);
  mklspace::mkl_free(Bcsr);
  mklspace::mkl_free(AJ);
  mklspace::mkl_free(BJ);
  mklspace::mkl_free(AI);
  mklspace::mkl_free(BI);
  debug("done\n%s","");

  return;
}

void check_status(mklspace::sparse_status_t status, char *label) {
  if (status != mklspace::SPARSE_STATUS_SUCCESS) {
    fprintf(stderr, "%s: unsuccessful: %d", label, status);
    fflush(stderr);
    exit(1);
  }
}

void mkl_spblas_spmmd(gsl_matrix *A, const gsl_matrix* B, gsl_matrix* C) {
  const MKL_INT an = A->size1;
  const MKL_INT am = A->size2;
  const MKL_INT bn = B->size1;
  const MKL_INT bm = B->size2;
  MKL_INT      ldc = C->size1;
  MKL_INT       NZ = 0;
  MKL_INT*      AI = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * (an+1),64);
  MKL_INT*      AE = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * (an),  64);
  MKL_INT*      AB = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * (an),  64);
  MKL_INT*      BI = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * (bn+1),64);
  MKL_INT*      BE = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * (bn),  64);
  MKL_INT*      BB = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * (bn),  64);
  MKL_INT     info;       // Output status from MKL

  // Note that this is not a general solver because NZ(A) is expected to equal NZ(B):
  double* Acsr;
  double* Bcsr;
  MKL_INT*  AJ;
  MKL_INT*  BJ;

#pragma omp parallel shared(NZ)
  {
#pragma omp for reduction(+:NZ) schedule(static)
    for (long i = 0; i < an*am; i++)
      if (A->data[i] != 0.0)
	NZ = NZ+1;
  }

  plainbug("a ");
  MKL_INT job[6] = {0, 0, 0, 2, NZ, 1};
  Acsr   = (double*)  mklspace::mkl_malloc(sizeof(double)  * NZ, 64);
  Bcsr   = (double*)  mklspace::mkl_malloc(sizeof(double)  * NZ, 64);
  AJ     = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * NZ, 64);
  BJ     = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * NZ, 64);

#pragma omp parallel num_threads(2)
  {
    #pragma omp task
    mklspace::mkl_ddnscsr(job, &an, &am, A->data, &am, Acsr, AJ, AI, &info);
    #pragma omp task
    mklspace::mkl_ddnscsr(job, &bn, &bm, B->data, &bm, Bcsr, BJ, BI, &info);
  }

  plainbug("b ");
  for (long i = 0; i < an; i++) {
    AB[i] = AI[i];
    AE[i] = AI[i+1];
  }
  for (long i = 0; i < bn; i++) {
    BB[i] = BI[i];
    BE[i] = BI[i+1];
  }

  plainbug("c ");
  mklspace::sparse_matrix_t Asparse;
  mklspace::sparse_matrix_t Bsparse;
  mklspace::sparse_status_t status;

#pragma omp parallel num_threads(2)
  {
    #pragma omp task
    status = mklspace::mkl_sparse_d_create_csr(&Asparse, mklspace::SPARSE_INDEX_BASE_ZERO, an, am, AB, AE, AJ, Acsr);
    #pragma omp task
    status = mklspace::mkl_sparse_d_create_csr(&Bsparse, mklspace::SPARSE_INDEX_BASE_ZERO, bn, bm, BB, BE, BJ, Bcsr);
  }

  check_status(status,"create_csr");
  const mklspace::sparse_operation_t op = mklspace::SPARSE_OPERATION_TRANSPOSE;
  struct mklspace::matrix_descr descr; //= { .type = SPARSE_MATRIX_TYPE_GENERAL };

  if (false) {
    plainbug("d(hinting) ");
    descr.type = mklspace::SPARSE_MATRIX_TYPE_GENERAL;
    //status = mklspace::mkl_sparse_set_memory_hint( Asparse, SPARSE_MEMORY_NONE );
    //check_status( status, "set_memory_hint a" );
    status = mklspace::mkl_sparse_set_mm_hint(Asparse, op, descr, mklspace::SPARSE_LAYOUT_ROW_MAJOR, an, 1);
    check_status( status, "set_mm_hint a" );
    status = mklspace::mkl_sparse_optimize( Asparse );
    check_status( status, "optimize a" );

    //status = mklspace::mkl_sparse_set_memory_hint( Bsparse, SPARSE_MEMORY_NONE );
    //check_status( status, "set_memory_hint b" );
    status = mklspace::mkl_sparse_set_mm_hint(Bsparse, op, descr, mklspace::SPARSE_LAYOUT_ROW_MAJOR, bn, 1);
    check_status( status, "set_mm_hint b" );
    status = mklspace::mkl_sparse_optimize( Bsparse );
    check_status( status, "optimize b" );
  }

  plainbug("e ");
  status = mklspace::mkl_sparse_d_spmmd(op, Asparse, Bsparse,
					mklspace::SPARSE_LAYOUT_ROW_MAJOR,
					C->data, ldc);
  check_status(status,"multiply(C := B^T * A)");

  plainbug("f ");
  mkl_sparse_destroy(Asparse);
  mkl_sparse_destroy(Bsparse);
  mklspace::mkl_free(Acsr);
  mklspace::mkl_free(Bcsr);
  mklspace::mkl_free(BI);
  mklspace::mkl_free(BJ);
  mklspace::mkl_free(BB);
  mklspace::mkl_free(BE);
  mklspace::mkl_free(AI);
  mklspace::mkl_free(AJ);
  mklspace::mkl_free(AB);
  mklspace::mkl_free(AE);

  plainbug(" done\n");
  return;
}



void compcol2csr31(mkl_spmat* S, gsl_spmatrix* C) {
  S->NZ = C->nz;
  S->n  = C->size1;
  S->m  = C->size2;
  S->Mi = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * (S->n )+1, 64);
  S->Mv = (double*)  mklspace::mkl_malloc(sizeof(double)  *  S->NZ,    64);
  S->Mj = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) *  S->NZ,    64);

  gsl_spmatrix* CT = gsl_spmatrix_alloc_nzmax(S->n,S->m,S->NZ,GSL_SPMATRIX_CCS);
  gsl_spmatrix_transpose_memcpy(CT,C);

  debug("S7%s\n","");
  for (long i = 0; i < S->NZ; i++) {
    S->Mv[i] = CT->data[i];
    S->Mj[i] = CT->i[i]+1;
  }
  debug("S6%s\n","");
  for (long i = 0; i < (S->n)+1; i++)
    S->Mi[i] = CT->p[i]+1;

  S->Mi[(S->n)+1] = S->NZ+1;
  
  gsl_spmatrix_free(CT);
 
  return;
}

int mkl_csradd(mkl_spmat* &A, mkl_spmat* B) {
  char    trans = 'N';
  MKL_INT req   = 1;
  double  beta  = 1.0;
  MKL_INT m     = A->m;
  MKL_INT n     = A->n;
  MKL_INT sort  = 0;
  MKL_INT nzmax = A->NZ + B->NZ;
  MKL_INT info;
  
  debug("A1%s\n","");
  mkl_spmat *C = (mkl_spmat*) malloc(sizeof(mkl_spmat));
  C->n  = n;
  C->m  = m;
  C->Mi = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * (n+1), 64);
  
  debug("A2%s\n","");
  mklspace::mkl_dcsradd(&trans, &req,  &sort, &m, &n, A->Mv, A->Mj, A->Mi,  &beta,
			B->Mv,  B->Mj, B->Mi, C->Mv,  C->Mj, C->Mi, &nzmax, &info);

  debug("A3%s\n","");
  C->NZ = C->Mi[n]-1;
  C->Mv = (double*)  mklspace::mkl_malloc(sizeof(double)  * C->NZ, 64);
  C->Mj = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * C->NZ, 64);

  req = 2;
  debug("A4%s\n","");
  mklspace::mkl_dcsradd(&trans, &req,  &sort, &m, &n, A->Mv, A->Mj, A->Mi,  &beta,
			B->Mv,  B->Mj, B->Mi, C->Mv,  C->Mj, C->Mi, &nzmax, &info);
  
  debug("A5%s\n","");
  /*
  // Copy result, C, into A:
  A->n  = C->n;
  A->m  = C->m;
  A->NZ = C->NZ;

  debug("A6%s\n","");
  for (long i = 0; i < n+1; i++)
    A->Mi[i] = C->Mi[i];

  debug("A7%s\n","");
  for (long i = 0; i < C->NZ; i++) {
    A->Mv[i] = C->Mv[i];
    A->Mj[i] = C->Mj[i];
  }
  */
  mkl_spmat_free(A);
  debug("A8%s\n","");

  A = C;
  debug("A9%s\n","");

  return info;
}

int mkl_csrmul(mkl_spmat* A, mkl_spmat* B, mkl_spmat* &C) {
  char   trans  = 'N';
  MKL_INT  req  = 1;
  MKL_INT n     = A->n;
  MKL_INT m     = A->m;
  MKL_INT k     = B->n;
  MKL_INT sort  = 7;
  MKL_INT nzmax = A->NZ + B->NZ;
  MKL_INT info;
  
  C     = (mkl_spmat*) malloc(sizeof(mkl_spmat));
  C->n  = n;
  C->m  = m;
  C->Mi = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * (m+1), 64);

  mklspace::mkl_dcsrmultcsr(&trans, &req, &sort, &m, &n, &k, A->Mv, A->Mj, A->Mi, 
			    B->Mv, B->Mj, B->Mi, C->Mv, C->Mj, C->Mi, &nzmax, &info);
  
  C->NZ = C->Mi[m]-1;
  C->Mv = (double*)  mklspace::mkl_malloc(sizeof(double)  * C->NZ, 64);
  C->Mj = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * C->NZ, 64);

  req=2;
  mklspace::mkl_dcsrmultcsr(&trans, &req, &sort, &m, &n, &k, A->Mv, A->Mj, A->Mi, 
			    B->Mv, B->Mj, B->Mi, C->Mv, C->Mj, C->Mi, &nzmax, &info);
  
  return info;
}

#endif

#ifdef SPARSE
void PrepareRegressionComponents(corpus_t*    corpus,    lda_seq*    seq,            unsigned int k,
    				 gsl_matrix*  W_phi,     gsl_matrix* W_phi_var,      gsl_vector*  d_phi_var_tmp,
                                 gsl_vector*  response,  gsl_vector* l_update_term2, gsl_matrix*  l_update_term1,
                                 gsl_matrix* W_phi_tmp, gsl_vector*  exp_h_tmp,
				 const size_t t, const int nmetafields, gsl_spmatrix*  tau, 
				 const bool write_regression, const char* root) {
#else
void PrepareRegressionComponents(corpus_t*    corpus,    lda_seq*     seq,            unsigned int k,
    				 gsl_matrix*  W_phi,     gsl_matrix*  W_phi_var,      gsl_vector*  d_phi_var_tmp,
                                 gsl_vector*  response,  gsl_vector*  l_update_term2, gsl_matrix*  l_update_term1,
				 gsl_matrix*  W_phi_tmp, gsl_vector*  exp_h_tmp,
				 const size_t t,   const int nmetafields,             gsl_matrix*  tau, 
				 const bool   write_regression,                 const char* root) {
#endif
  debug("DEBUG 4%s\n","");
  gsl_blas_dgemv(CblasTrans, 1.0, W_phi, response, 0.0, l_update_term2);
  debug("DEBUG 5%s\n","");

  // Set up the transformation matrix.
  // First, set up W_phi^T \Lambda W_phi.
  gsl_matrix_memcpy(W_phi_tmp, W_phi);

#pragma omp parallel 
{
  gsl_vector_view col;
#pragma omp for schedule(static)
  for (int d=0; d < corpus->ndocs; ++d) {
    col = gsl_matrix_column(W_phi_tmp, d);
    gsl_vector_mul(&col.vector, exp_h_tmp);
  }
  }

  debug("DEBUG 6%s\n","");

#ifdef MKL
#ifdef TESTDGEMM
  if (t==1 || t==20 || t==40 || t==60 || t==80) {
    if (k < 5) {
      double startt = omp_get_wtime();
      mkl_cblas_dgemm(CblasTrans, CblasNoTrans, 1.0, W_phi_tmp, W_phi, 0.0, l_update_term1);
      debug("DEBUG 7 **** [SPEC] MKL dgemm took %f\n", omp_get_wtime()-startt);
    }
    else if (k >= 5 && k < 10) {
      double startt = omp_get_wtime();
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W_phi_tmp, W_phi, 0.0, l_update_term1);
      debug("DEBUG 7 **** [SPEC] GSL dgemm took %f\n", omp_get_wtime()-startt);
    }
    else if (k >= 10 && k < 15) {
      double startt = omp_get_wtime();
      mkl_spblas_dgemm(W_phi_tmp, W_phi, l_update_term1);
      debug("DEBUG 7 **** [SPEC] Sparse MKL dgemm took %f\n", omp_get_wtime()-startt);
    }
    else {
      double startt = omp_get_wtime();
      mkl_spblas_spmmd(W_phi_tmp, W_phi, l_update_term1);
      debug("DEBUG 7 **** [SPEC] Sparse MKL spmmd took %f\n", omp_get_wtime()-startt);
    }
  }
  else {
    double startt = omp_get_wtime();
    mkl_spblas_spmmd(W_phi_tmp, W_phi, l_update_term1);
    debug("DEBUG 7 **** Sparse MKL spmmd took %f\n", omp_get_wtime()-startt);
  }
#else
  double startt = omp_get_wtime();
  mkl_spblas_spmmd(W_phi_tmp, W_phi, l_update_term1);
  debug("DEBUG 7 **** Sparse MKL spmmd took %f\n", omp_get_wtime()-startt);
#endif

#else
  double startt = omp_get_wtime();
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W_phi_tmp, W_phi, 0.0, l_update_term1);
  debug("DEBUG 7 **** GSL dgemm took %f\n", omp_get_wtime()-startt);
#endif
  
  // d_phi_var_tmp = W_phi_var^T * exp_h_tmp
  gsl_blas_dgemv(CblasTrans, 1.0, W_phi_var, exp_h_tmp, 0.0, d_phi_var_tmp);
  
  debug("DEBUG 8%s\n","");

  // rDIM update for \ell. Equation 21 in YH's derivation
  // Added by AG
  if (FLAGS_model == "rdim" && tau != NULL) {
    gsl_vector      *extrinsic_regression_term = gsl_vector_alloc(corpus->ndocs);
    gsl_vector_view  mu_k                      = gsl_matrix_column(seq->mu, k);

    // extrinsic_regression_term = mu_k * tau_t
#ifdef SPARSE
    debug("DEBUG 9%s\n","");
    gsl_spblas_dgemv(CblasNoTrans, 1.0, tau, &mu_k.vector, 0.0, extrinsic_regression_term);
#else
    gsl_blas_dgemv(CblasNoTrans, 1.0, tau, &mu_k.vector, 0.0, extrinsic_regression_term);
#endif
    
    debug("DEBUG 10%s\n","");
    // Add extrinsic_regression_term to l_update_term2:
    gsl_vector_add(l_update_term2, extrinsic_regression_term);
    
    // If this is the final iteration, we want to write
    // the \mu_k\tau_t to a file
    if (write_regression) {
      char name[400];
      sprintf(name, "%s/inf_update-regression_vector-%d.dat", root, (int) t);
      FILE* f = fopen(name, "w");
      gsl_vector_fprintf(f, extrinsic_regression_term, "%f");
      fclose(f);
    }
    
    debug("DEBUG 11%s\n","");
    gsl_vector_free(extrinsic_regression_term);
  } // End rDIM update
  
  // Next, add elements to the diagonal of l_update_term1.
#pragma omp parallel for schedule(static) num_threads(FLAGS_threads)
  for (unsigned int d=0; d < corpus->ndocs; ++d) {
    if (FLAGS_debug){
	if (d==0) {
	   debug("We have %i OMP threads\n",omp_get_num_threads());
	}
    }
    double value = gsl_matrix_get(l_update_term1, d, d);
    value += (seq->topic[k]->chain_variance / (FLAGS_sigma_d * FLAGS_sigma_d));

    // sgerrish: Is this supposed to be multiplied by anything?
    value += gsl_vector_get(d_phi_var_tmp, d);
    gsl_matrix_set(l_update_term1, d, d, value); // Note, only setting diagonal
  }
  debug("DEBUG 12%s\n","");
}

void SetExpHTmp(lda_seq*     seq, const corpus_seq_t* data,      unsigned int t,      unsigned int k,
	        gsl_vector*  exp_h_tmp, gsl_vector*  zw_tmp, gsl_vector** exp_i_tmp) {
  gsl_vector_set_zero(exp_h_tmp);

  for (int i = t; i < seq->nseq; ++i) {
    gsl_vector_view mean_i_current = gsl_matrix_column(seq->topic[k]->e_log_prob, i);
    gsl_vector_view var_i_current  = gsl_matrix_column(seq->topic[k]->variance, i + 1);

    // Set up exp_h_tmp.
    gsl_vector_memcpy(zw_tmp, &var_i_current.vector);
    gsl_vector_sub(zw_tmp, &mean_i_current.vector);
    gsl_vector_scale(zw_tmp, 2.0);

    // Set up exp_i_tmp.
    gsl_vector_memcpy(exp_i_tmp[i], &var_i_current.vector);
    gsl_vector_scale(exp_i_tmp[i], 0.5);
    gsl_vector_sub(exp_i_tmp[i], &mean_i_current.vector);

#pragma omp parallel for schedule(static) num_threads(FLAGS_threads)
    for (int n=0; n < data->nterms; ++n) {
      gsl_vector_set(zw_tmp, n, exp(gsl_vector_get(zw_tmp, n)));
      gsl_vector_set(exp_i_tmp[i], n, exp(gsl_vector_get(exp_i_tmp[i], n)));
    }
    
    gsl_vector_scale(zw_tmp, scaled_influence[i - t] * scaled_influence[i - t]);
    gsl_vector_add(exp_h_tmp, zw_tmp);
  }
}

// Updater for the influence variable in DIM
double update_inf_var_fixed(lda_seq* seq, const corpus_seq_t* data, gsl_matrix**   phi, 
			    size_t     t, const char*         root, int dump_doc_stats,
			    const bool write_regression,     const bool        heldout) {
  double lhood = 0.0;

  // Note that we're missing a suspicious factor of -2 on the document
  // weights currently.  We won't worry about that for now (since
  // we're still experimenting), but it should soon be fixed. // SG
  corpus_t*  corpus = data->corpus[t];

  if (!corpus->ndocs) 
    return lhood;
  
  inf_var* influence = seq->influence; // Comes from training, never test...
  if (t != TIME && 0) // Whaaa!?? // AG
    return lhood;

  // We need access to the following:
  gsl_matrix* documents_topics_t              = NULL;
  gsl_matrix* renormalized_documents_topics_t = NULL;
  
  if (!heldout) {
    documents_topics_t              = influence->doc_weights[t]; 
    renormalized_documents_topics_t = influence->renormalized_doc_weights[t];
  } else {
    documents_topics_t              = gsl_matrix_calloc(corpus->ndocs, seq->ntopics);
    renormalized_documents_topics_t = gsl_matrix_calloc(corpus->ndocs, seq->ntopics);
  }

  gsl_matrix
    *W_phi          = gsl_matrix_calloc(seq->nterms,   corpus->ndocs),
    *W2_phi2        = gsl_matrix_calloc(seq->nterms,   corpus->ndocs),
    *W_phi_var      = gsl_matrix_calloc(seq->nterms,   corpus->ndocs),
    *W_phi_tmp      = gsl_matrix_alloc(seq->nterms,   corpus->ndocs),
    *l_update_term1 = gsl_matrix_alloc(corpus->ndocs, corpus->ndocs);
  gsl_vector
    *xd_tmp         = gsl_vector_alloc(corpus->ndocs),
    *xd_tmp2        = gsl_vector_alloc(corpus->ndocs),
    *yw_tmp         = gsl_vector_alloc(seq->nterms),
    *zw_tmp         = gsl_vector_alloc(seq->nterms),
    *terms_inc_tmp  = gsl_vector_alloc(seq->nterms);
  gsl_vector
    **exp_i_tmp     = new gsl_vector*[seq->nseq];

  for (int i=0; i < seq->nseq; ++i)
    exp_i_tmp[i] = gsl_vector_alloc(seq->nterms);
  
  gsl_vector
    *response       = gsl_vector_alloc(seq->nterms), // Get's zero'd every iteration
    *exp_h_tmp      = gsl_vector_alloc(seq->nterms),
    *l_update_term2 = gsl_vector_alloc(corpus->ndocs),
    *d_phi_var_tmp  = gsl_vector_alloc(corpus->ndocs);
  double
    *total_terms            = (double*) malloc(sizeof(double) * corpus->ndocs),
    *renormalization_totals = (double*) malloc(sizeof(double) * corpus->ndocs);

  // Line 2 in rDIM algorithm

  // MPI NOTE: This loop is is what I think could be distributed by rank:
  // There is OpenMP threading beneatch this
  for (size_t k=0; k < documents_topics_t->size2; ++k) { // For each topic
    // Set up W_phi and W_phi_var.
    // Line 3 in rDIM algorithm
    debug("DEBUG a0%s\n","");

#pragma omp parallel for schedule(static) num_threads(FLAGS_threads)
    for (int d=0; d < corpus->ndocs; ++d) { // For every document
      doc_t* doc = corpus->doc[d];
      total_terms[d] = 0.0;

      for (int n=0; n < doc->nterms; ++n)
	total_terms[d]+=doc->count[n];

      // Assumed FLAGS_normalize_docs == "normalize" || FLAGS_normalize_docs == "occurrence"):
      renormalization_totals[d] = 1.0;

      for (int n=0; n < doc->nterms; ++n) { // For each word in the vocab of d
	// Phi, for the doc's term n and topic k.

	// Reordered some algebra in here -- might break unit-tests for associative math DPOPS:

	double  // Tried to consolidate these DPOs, but they aren't IEEE-transitive
	  phi_d_n_k      = gsl_matrix_get(phi[d], n, k), // Actually /phi_t_d_n_k
	  number_terms   = ((double) doc->count[n] / (double) total_terms[d]),
	  W_phi_term     = number_terms * phi_d_n_k,
	  phi_d_n_k2     = phi_d_n_k * phi_d_n_k,
	  number_terms2  = number_terms * number_terms;
	double 
	  W2_phi2_term   = number_terms2 *  phi_d_n_k2,
	  W_phi_var_term = number_terms2 * (phi_d_n_k - phi_d_n_k2);
	
	  gsl_matrix_set(W_phi, doc->word[n],     d, W_phi_term);
	  gsl_matrix_set(W2_phi2, doc->word[n],   d, W2_phi2_term);
	  gsl_matrix_set(W_phi_var, doc->word[n], d, W_phi_var_term); 
      }
    }
    
    gsl_vector_view renormalized_document_weights = gsl_matrix_column(renormalized_documents_topics_t, k);    
    
    // Now, with w_phi_var, etc. set, determine \sum_{i=t}^{T-1} r(...) h(t, i)
    SetExpHTmp(seq, data, t, k, exp_h_tmp, zw_tmp, exp_i_tmp);
    debug("DEBUG a3%s\n","");
    
    // Next, set up the weighted response,
    // \exp_{-m + v / 2) \circ (m_{t+1} - m_t + v_t}).
    // Here we also subtract the current l's contribution to influence_sum_lgl.
    gsl_vector_view w_phi_l_t = gsl_matrix_column(seq->topic[k]->w_phi_l, t);
    gsl_vector_set_zero(response);
    
    gsl_vector_view total_influence_time_i;
    gsl_vector_view mean_i_current;
    gsl_vector_view mean_i_next;
    gsl_vector_view var_i_current;
    
    debug("DEBUG a4%s\n","");
      
    // Tried hard to OpenMP this, but couldn't get it to replicate FPOs...
    for (int i = t; i < seq->nseq - 1; ++i) { // For every future time-slice
      total_influence_time_i = gsl_matrix_column(seq->influence_sum_lgl[k], i);
      gsl_vector_memcpy(zw_tmp, &w_phi_l_t.vector);
      gsl_vector_scale( zw_tmp,  scaled_influence[i - t]);
	
      gsl_vector_sub(&total_influence_time_i.vector, zw_tmp);
	
      gsl_vector_memcpy(zw_tmp, &total_influence_time_i.vector);
      gsl_vector_mul(zw_tmp,  exp_i_tmp[i]);
	
      mean_i_current = gsl_matrix_column(seq->topic[k]->e_log_prob, i);
      mean_i_next    = gsl_matrix_column(seq->topic[k]->e_log_prob, i + 1);
      var_i_current  = gsl_matrix_column(seq->topic[k]->variance,   i + 1);
	
      gsl_vector_memcpy(terms_inc_tmp, &mean_i_next.vector);
      gsl_vector_sub(   terms_inc_tmp, &mean_i_current.vector);
      gsl_vector_add(   terms_inc_tmp, &var_i_current.vector);
      gsl_vector_sub(   terms_inc_tmp,  zw_tmp);
      gsl_vector_mul(   terms_inc_tmp,  exp_i_tmp[i]);
	
      gsl_vector_scale( terms_inc_tmp,  scaled_influence[i - t]);
	
      gsl_vector_add(response, terms_inc_tmp);
    }

    debug("DEBUG a5%s\n","");
    // If were in rDIM, this function will include the mu_k * tau_t component:
    // MPI NOTE: This is a bear:
    PrepareRegressionComponents(corpus, seq, k, W_phi, W_phi_var, d_phi_var_tmp, 
				response, l_update_term2, l_update_term1, 
			        W_phi_tmp, exp_h_tmp, t, data->nmetafields, 
				(FLAGS_model == "rdim") ? data->tau[t] : NULL,
				write_regression, root);
    
    // Finally, solve for the document weights d!
    // Solves the linear system Ax = b for x:
    // In rDIM, l_update_term2 will have the \mu_k\tau_t component added from above // AG

    gsl_vector *document_weights = gsl_vector_alloc(documents_topics_t->size1);

    // This is overwritten by the solve below. Why was this here: //AG:
    //gsl_matrix_get_col(document_weights, documents_topics_t, k);

    debug("DEBUG 14%s\n","");

    double startt = omp_get_wtime();
    // Note: l_update_term1 is not symmetric and indefinite (pardiso mtype=11)
    //MKLSolve_lapacke_single(l_update_term1, l_update_term2, document_weights); // Preserves l_update_term1
    MKLSolve_lapacke_comp(l_update_term1, l_update_term2, document_weights);     // Clobbers  l_update_term1
    //LUSolve(l_update_term1, l_update_term2, document_weights);                 // Preserves l_update_term1
    //CholeskySolve(l_update_term1, l_update_term2, document_weights);           // Clobbers  l_update_term1
    debug("MKL solve took %f\n", omp_get_wtime()-startt);

    debug("DEBUG 15%s\n","");

    if (FLAGS_save_time == -1 || FLAGS_save_time == t) { // && t > 0) {
      outlog("Saving data at t = %d.\n", (int) t);
      
      char name[400];
      sprintf(name, "%s%ld_%ld_weighted_document_terms.dat", root, k, t);
      FILE* f = fopen(name, "w");
      params_write_sparse_gsl_matrix(f, (char* const) "W_phi", W_phi);
      fclose(f);
      
      sprintf(name, "%s%ld_%ld_weighted_document_terms_var.dat", root, k, t);
      f = fopen(name, "w");
      params_write_sparse_gsl_matrix(f, (char* const) "W_phi_var", W_phi_var);
      fclose(f);
      
      sprintf(name, "%s%ld_%ld_phi.dat", root, k, t);
      f = fopen(name, "w");
      params_write_gsl_matrix(f, (char* const) "phi 0", phi[0]);
      fclose(f);
      
      sprintf(name, "%s%ld_%ld_weighted_document_terms_sq.dat", root, k, t);
      f = fopen(name, "w");
      params_write_sparse_gsl_matrix(f, (char* const) "W2_phi2", W2_phi2);
      fclose(f);
      
      /* sprintf(name, "%s%ld_%ld_weighted_response.dat", root, k, t); */
      /* f = fopen(name, "w"); */
      /* params_write_gsl_vector_multiline(f, (char* const) "weighted_response", l_update_term2); */
      /* fclose(f); */
      
      sprintf(name, "%s%ld_%ld_response.dat", root, k, t);
      f = fopen(name, "w");
      params_write_gsl_vector_multiline(f, (char* const) "response", response);
      fclose(f);
      
      sprintf(name, "%s%ld_%ld_exp_h.dat", root, k, t);
      f = fopen(name, "w");
      params_write_gsl_vector_multiline(f, (char* const) "exp_h", exp_h_tmp);
      fclose(f);
            
      DumpTimeDocTopicStats(root, t, corpus, phi);
    }
    /* We really don't need this mid-iteration: */
    /* else if (FLAGS_save_time == -1) // First, dump phi's for the top topics. */
    /*   DumpTimeDocTopicStats(root, t, corpus, phi); */

    debug("DEBUG 17%s\n","");
    outlog("Done updating topic %ld, time %ld.\n", k, t);

    // AG: ?Applying the normalizations?
#pragma omp parallel for schedule(static) num_threads(FLAGS_threads)
    for (int d = 0; d < document_weights->size; ++d) 
      gsl_vector_set(&renormalized_document_weights.vector, d,
		     vget(document_weights, d) * renormalization_totals[d]);

    // Now copy this and several products to the sslm_var object.
    gsl_vector_view w_phi_l = gsl_matrix_column(seq->topic[k]->w_phi_l, t);

    // Copies the product of document_weights and w_phi_l into W_phi
    gsl_blas_dgemv(CblasNoTrans, 1.0, W_phi, document_weights, 0.0, &w_phi_l.vector);
    debug("DEBUG 18%s\n","");
    
    // Copy this value back into lgl.
   
#pragma omp parallel num_threads(FLAGS_threads)
    {
      gsl_vector* zw_tmp_g = gsl_vector_alloc(seq->nterms);
#pragma omp for schedule(static)
    for (int i=t; i < seq->nseq - 1; ++i) {
      if (i==t) {
        debug("OMP DEBUG G (%i threads)\n",omp_get_num_threads());
      }
      gsl_vector_view total_influence_time_i = gsl_matrix_column(seq->influence_sum_lgl[k], i);
      
      gsl_vector_memcpy(zw_tmp_g, &w_phi_l.vector);
      gsl_vector_scale(zw_tmp_g, scaled_influence[i - t]);
#pragma omp critical
      {
	gsl_vector_add(&total_influence_time_i.vector, zw_tmp_g);
      }
    }
    gsl_vector_free(zw_tmp_g);
    } // GOMP

    debug("DEBUG 19%s\n","");
    // Keep track of the term we need to add to m_update_coeff.
    gsl_vector_memcpy(terms_inc_tmp, &w_phi_l.vector);
    gsl_vector_mul(terms_inc_tmp, &w_phi_l.vector);

    // Copy and square the document weights vector.
#pragma omp parallel for schedule(static) num_threads(FLAGS_threads)
    for (int i = 0; i < xd_tmp->size; ++i) {
      double value = gsl_vector_get(document_weights, i);
      value = value * value + FLAGS_sigma_l * FLAGS_sigma_l;
      gsl_vector_set(xd_tmp, i, value);
      gsl_vector_set(xd_tmp2, i, FLAGS_sigma_l * FLAGS_sigma_l);
    }

    debug("DEBUG 20%s\n","");
    gsl_blas_dgemv(CblasNoTrans, 1.0, W_phi_var, xd_tmp, 0.0, yw_tmp);
    gsl_vector_add(terms_inc_tmp, yw_tmp);

    debug("DEBUG 21%s\n","");
    
    gsl_blas_dgemv(CblasNoTrans, 1.0, W2_phi2, xd_tmp2, 0.0, yw_tmp);
    gsl_vector_add(terms_inc_tmp, yw_tmp);
    debug("DEBUG 22%s\n","");

    // Store an update coefficient for the beta updates. (not the rDIM coefficients)
#pragma omp parallel num_threads(FLAGS_threads)
    {
      gsl_vector* yw_tmp_g = gsl_vector_alloc(seq->nterms);

#pragma omp for schedule(static)
    for (int i = t; i < seq->nseq; ++i) {
      gsl_vector_view m_update_coeff_h = gsl_matrix_column(seq->topic[k]->m_update_coeff, i);
      gsl_vector_view m_update_coeff_g = gsl_matrix_column(seq->topic[k]->m_update_coeff_g, i);
      if (t == 0) {
	gsl_vector_set_zero(&m_update_coeff_h.vector);
	gsl_vector_set_zero(&m_update_coeff_g.vector);       
      }

      gsl_vector_memcpy(yw_tmp_g, terms_inc_tmp);
      gsl_vector_scale(yw_tmp_g, scaled_influence[i-t]);

      gsl_vector_add(&m_update_coeff_h.vector, yw_tmp_g);

      gsl_vector_memcpy(yw_tmp_g, &w_phi_l.vector);
      gsl_vector_scale(yw_tmp_g, scaled_influence[i - t]);

      gsl_vector_add(&m_update_coeff_g.vector, yw_tmp_g);
    }
    gsl_vector_free(yw_tmp_g);
    } // GOMP

    debug("DEBUG 24%s\n","");
#pragma omp parallel shared(lhood) num_threads(FLAGS_threads)
    {
#pragma omp for reduction(+:lhood) schedule(static)
      for (int i = 0; i < corpus->ndocs; ++i) {
	double value = gsl_vector_get(document_weights, i);
	value = (-(value * value + FLAGS_sigma_l * FLAGS_sigma_l)
		 / (2.0 * FLAGS_sigma_d * FLAGS_sigma_d)
		 - 0.5 * log(2 * PI) - log(FLAGS_sigma_d * FLAGS_sigma_d));
	lhood = lhood + value;
      }
    }
  
    gsl_matrix_set_col(documents_topics_t, k, document_weights);
    gsl_vector_free(document_weights);  

    debug("DEBUG 25%s\n","");
  }

  free(total_terms);
  free(renormalization_totals);

  gsl_matrix_free(W_phi);
  gsl_matrix_free(W_phi_tmp);
  gsl_matrix_free(W2_phi2);
  gsl_matrix_free(W_phi_var);
  gsl_matrix_free(l_update_term1);
  gsl_vector_free(exp_h_tmp);
  gsl_vector_free(response);
  gsl_vector_free(terms_inc_tmp);

  for (int i=0; i < seq->nseq; ++i)
    gsl_vector_free(exp_i_tmp[i]);
  
  delete[] exp_i_tmp;
  gsl_vector_free(l_update_term2);
  gsl_vector_free(d_phi_var_tmp);

  gsl_vector_free(xd_tmp);
  gsl_vector_free(xd_tmp2);
  gsl_vector_free(yw_tmp);
  gsl_vector_free(zw_tmp);

  if (heldout) {
    gsl_matrix_free(documents_topics_t);
    gsl_matrix_free(renormalized_documents_topics_t);
  }
  
  return lhood;
}

void make_lda_from_seq_slice(lda* lda_m, lda_seq* lda_seq_m, int time) {
     // set lda model topics
     // !!! note: we should be able to point to the view...

     for (int k = 0; k < lda_seq_m->ntopics; k++) {
       // get topic
       gsl_vector s = gsl_matrix_column(lda_seq_m->topic[k]->e_log_prob, time).vector;
       gsl_vector d = gsl_matrix_column(lda_m->topics, k).vector;
       gsl_blas_dcopy(&s, &d);
     }

     gsl_blas_dcopy(lda_seq_m->alpha, lda_m->alpha);
}

// This is a consolidation of g{1-5}_alloc()s in the DTM code
static void consolidated_g_alloc(gsl_matrix** g, gsl_matrix** g3_matrix, 
				 gsl_matrix** g4_matrix, gsl_matrix** g5_matrix, 
				 lda_seq* model, const corpus_seq_t* data, int time) {

  *g         = gsl_matrix_alloc(model->nterms, model->ntopics);
  *g3_matrix = gsl_matrix_alloc(model->nterms, model->ntopics);
  *g4_matrix = gsl_matrix_alloc(model->nterms, model->ntopics);
  *g5_matrix = gsl_matrix_alloc(model->nterms, model->ntopics);  

  debug("FLAGS_kthreads=%i,FLAGS_threads=%i\n",FLAGS_kthreads,FLAGS_threads);

#pragma omp parallel num_threads(FLAGS_kthreads)
  {

  double 
    exp_m, m, m_next, variance_first,
    variance, total2, total3, total4, exp_m_scaled, 
    w_phi_l, iot;

//#pragma omp for schedule(dynamic) // Old way was dynamic
#pragma omp for schedule(static)
  for (int k = 0; k < model->ntopics; ++k) { // 10~30
    for (int w = 0; w < model->nterms; ++w) { // ~20,000
      variance_first = mget(model->topic[k]->variance,   w, time);
      m              = mget(model->topic[k]->e_log_prob, w, time);
      exp_m          = exp(-m + variance_first / 2.0);

      gsl_matrix_set(*g, w, k, (scaled_influence[0] * -variance_first * exp_m));
      
      total2 = 0.0;
      total3 = 0.0;
      total4 = 0.0;

      for (int i=time; i < model->nseq - 1; ++i) { // ~90
	iot = 0.0;
	
	for (int j = 0; j < i; ++j) {
	  exp_m    = exp(-mget(model->topic[k]->e_log_prob, w, j) + mget(model->topic[k]->variance, w, j) / 2.0);
	  iot     += (mget(model->topic[k]->w_phi_l, w, j) * scaled_influence[i - j] * exp_m);

	  w_phi_l  = mget(model->topic[k]->w_phi_l, w, j);
	  total3  += exp_m_scaled * w_phi_l * scaled_influence[i - j];
	}
	
	// This substraction is probably faster than testing j!=time above:
	iot -= (mget(model->topic[k]->w_phi_l, w, time) * scaled_influence[i - time] * exp_m);
	
	m       = mget(model->topic[k]->e_log_prob, w, i);
	m_next  = mget(model->topic[k]->e_log_prob, w, i + 1);
	gsl_matrix_set(*g, w, k, mget(*g, w, k) + (scaled_influence[i - time] * 
						   (m_next - m - iot)));
	
	variance      = mget(model->topic[k]->variance, w, i + 1);
	exp_m         = exp(-m + variance / 2.0);
	total2       += (scaled_influence[i - time] * exp_m * (m_next - m + variance));
	
	exp_m         = exp(-2.0 * m + 2.0 * variance);
	exp_m_scaled  = exp_m *  scaled_influence[i - time];
	total4       += exp_m * (scaled_influence[i - time] * scaled_influence[i - time]);
      }
      
      exp_m = exp(-m + mget(model->topic[k]->variance, w, time) / 2.0);
      
      gsl_matrix_set(*g,         w, k, mget(*g, w, k) * exp_m);
      
      gsl_matrix_set(*g3_matrix, w, k, total2);
      gsl_matrix_set(*g4_matrix, w, k, total3);
      gsl_matrix_set(*g5_matrix, w, k, total4);
    }
  }
  } // GOMP?

  return;
}

 /*
  * compute the likelihood of a sequential corpus under an LDA seq
  * model. return the likelihood bound.
  *
  */
static void InferDTMSeq(const int K, unsigned int iter, unsigned int last_iter, 
			const corpus_seq_t* data, gsl_matrix* gammas,
			gsl_matrix* lhoods, lda* lda_model, lda_post* post, 
			lda_seq* model, gsl_matrix** suffstats, double* bound) {
  int doc_index = 0;

  for (int t = 0; t < data->len; t++) {
    // Prepare coefficients for the phi updates. This change is relatively painless.
    make_lda_from_seq_slice(lda_model, model, t);
    int ndocs = data->corpus[t]->ndocs;

    for (int d = 0; d < ndocs; d++) {
      gsl_vector gam   = gsl_matrix_row(gammas, doc_index).vector;
      gsl_vector lhood = gsl_matrix_row(lhoods, doc_index).vector;
      post->gamma = &gam;
      post->doc   = data->corpus[t]->doc[d];
      post->lhood = &lhood;
      double doc_lhood;

      // For now, only do the standard, phi-based update.
      if (iter == 0)
	doc_lhood = fit_lda_post_unrolled(iter, d, t, post, NULL,  NULL, NULL, NULL, NULL);
      else
	doc_lhood = fit_lda_post_unrolled(iter, d, t, post, model, NULL, NULL, NULL, NULL);

      if (suffstats != NULL)
	update_lda_seq_ss(t, data->corpus[t]->doc[d], post, suffstats);

      *bound += doc_lhood;
      doc_index++;
    }
  }
}

// Main E-step for DIM (as oppoed to DTM)
static void InferDIMSeq(const int   K,         unsigned int        iter,      unsigned int last_iter, 
			const char* file_root, const corpus_seq_t* data,      gsl_matrix*  gammas,    
			gsl_matrix* lhoods,    lda*                lda_model, lda_post*    post,
			lda_seq*    model,     gsl_matrix**        suffstats, double*      bound,
			const bool  heldout) {
  int doc_index = 0;

  for (int t = 0; t < data->len; t++) { // For every time-slice, t
    outlog("DIM E step for iteration=%d, time=%d\n", iter, t);
    omp_set_max_active_levels(2);
    omp_set_nested(1);
    //mklspace::mkl_set_dynamic(0); // Old way had mkl dyn=0 here
    omp_set_dynamic(0);

    // Prepare coefficients for the phi updates. This change is relatively painless.
    gsl_matrix *g, *g3_matrix, *g4_matrix, *g5_matrix;
    consolidated_g_alloc(&g, &g3_matrix, &g4_matrix, &g5_matrix, model, data, t);
    make_lda_from_seq_slice(lda_model, model, t);

    int          ndocs = data->corpus[t]->ndocs;
    gsl_matrix** phi_t = (gsl_matrix**) malloc(ndocs * sizeof(gsl_matrix*));
    double     doc_lhood;
    
    debug("DEBUG D%s\n","");
    double timed = omp_get_wtime();
    for (int d = 0; d < ndocs; d++) { // For every document in time t
      gsl_vector gam   =  gsl_matrix_row(gammas, doc_index).vector;
      gsl_vector lhood =  gsl_matrix_row(lhoods, doc_index).vector;
      post->gamma      = &gam;
      post->doc        =  data->corpus[t]->doc[d];
      post->lhood      = &lhood;

      if (iter == 0)
	doc_lhood = fit_lda_post_unrolled(iter, d, t, post, NULL, NULL, NULL, NULL, NULL);
      else
	doc_lhood = fit_lda_post_unrolled(iter, d, t, post, model, g, g3_matrix, g4_matrix, g5_matrix);

      if (suffstats != NULL) 
	update_lda_seq_ss(t, data->corpus[t]->doc[d], post, suffstats);
      
      phi_t[d]                 = gsl_matrix_alloc(post->doc->nterms, K);
      gsl_matrix_view phi_view = gsl_matrix_submatrix( post->phi, 0, 0, post->doc->nterms, K);

      gsl_matrix_memcpy(phi_t[d], &phi_view.matrix);

      *bound += doc_lhood;
      doc_index++;
    }

    debug("DEBUG E0%s\n","");
    debug("DEBUG E (D took %f)\n", omp_get_wtime()-timed);

    mklspace::mkl_thread_free_buffers(); // Old way didn't have this...
    omp_set_max_active_levels(2);
    omp_set_nested(1);
    //mklspace::mkl_set_dynamic(0); // Old way had mkl dyn=0 and omp dyn=1 here
    //omp_set_dynamic(1);

    if (t < data->len - 1) { // If this is not the last time-slice
      if (FLAGS_model == "fixed" || FLAGS_model == "rdim") {
        double l_bound = update_inf_var_fixed(model, data, phi_t, t, file_root,
					      last_iter || iter >= FLAGS_lda_sequence_min_iter,
					      (iter == last_iter), heldout);
	*bound += l_bound;
      }
    } else { // added by YH
      corpus_t* corpus = data->corpus[t]; // added by YH
      DumpTimeDocTopicStats(file_root, t, corpus, phi_t); // added by YH
    }
    
    for (int d=0; d < ndocs; ++d) 
      gsl_matrix_free(phi_t[d]);

    free(phi_t);

    gsl_matrix_free(g);
    gsl_matrix_free(g3_matrix);
    gsl_matrix_free(g4_matrix);
    gsl_matrix_free(g5_matrix);
  }
  
  outlog("%s\n", "E-Step complete");
}

// AG: This is the bulk of the E-step
double lda_seq_infer(lda_seq*    model,  const corpus_seq_t* data, gsl_matrix** suffstats, gsl_matrix* gammas,
		     gsl_matrix* lhoods, unsigned int        iter, unsigned int last_iter, const char* file_root,
		     const bool  heldout) {
  int K          = model->ntopics, W = model->nterms;
  outlog("In lda_seq_infer(): K=%d W=%d data->nterms=%d, data->max_nterms=%d time-slices=%d\n", 
	 K, W, data->nterms, data->max_nterms,data->len);
  double bound   = 0.0;
  lda* lda_model = new_lda_model(K, W);
  lda_post post;

  post.phi     = gsl_matrix_calloc(data->max_nterms, K);
  post.log_phi = gsl_matrix_calloc(data->max_nterms, K);
  post.model   = lda_model;
  
  if (FLAGS_model == "fixed" || FLAGS_model == "rdim") { // DIM
    InfluenceTotalFixed(model); 
    InferDIMSeq(K,      iter,       last_iter, file_root, data,       gammas, 
		lhoods, lda_model, &post,      model,     suffstats, &bound, heldout); 
  } else if (FLAGS_model == "dtm") { // DTM
    InferDTMSeq(K,      iter,       last_iter, data,  gammas, 
		lhoods, lda_model, &post,      model, suffstats, &bound);
  } else {
    outlog("Error. Unknown model '%s'.\n", FLAGS_model.c_str());
    exit(1);
  }

  gsl_matrix_free(post.phi);
  gsl_matrix_free(post.log_phi);
  free_lda_model(lda_model);

  return(bound);
}

/*
 * fit an lda sequence model: (ie. DTM)
 *
 * . for each time period
 * .     set up lda model with E[log p(w|z)] and \alpha
 * .     for each document
 * .         perform posterior inference
 * .         update sufficient statistics/likelihood
 * . maximize topics
 */
double fit_lda_seq(           lda_seq* model,   const corpus_seq_t* data, 
		   const corpus_seq_t* heldout, const         char* file_root) {
  int K = model->ntopics, W = model->nterms, k;

  if (FLAGS_model == "rdim") {
    outlog("Because this is the first iteration, making Eq 32 term1%s\n","");
    //    make_term1_csr3(data, model);
    make_term1(data, model);
  }
  
  gsl_matrix* topic_suffstats[K];
  gsl_matrix* gammas;
  gsl_matrix* lhoods;
  gsl_matrix* heldout_gammas = NULL;
  gsl_matrix* heldout_lhoods = NULL;
  
  double 
    bound,
    heldout_bound,
    convergence,
    old_bound;
  
  int          iter;
  short        final_iters_flag;
  unsigned int last_iter;
  
  FILE*        em_log;
  char name[400];
  char root[400];
  char em_log_filename[400];
  
  // If we're not checkpointing, go ahead an initialize:
  if (strcmp("", FLAGS_checkpoint_recover.c_str()) == 0) {
    // initialize sufficient statistics
    for (k = 0; k < K; k++)
      topic_suffstats[k] = gsl_matrix_calloc(W, data->len);

    // set up variables
    gammas = gsl_matrix_calloc(data->ndocs, K);
    lhoods = gsl_matrix_calloc(data->ndocs, K+1); // AG: why K+1?

    if (heldout != NULL) {
      heldout_gammas = gsl_matrix_calloc(heldout->ndocs, K);
      heldout_lhoods = gsl_matrix_calloc(heldout->ndocs, K+1); // Why +1 here?
    }

    bound         = 0;
    heldout_bound = 0;
    convergence   = LDA_SEQ_EM_THRESH + 1;

    // run EM
    iter             = 0;
    final_iters_flag = 0;
    last_iter        = 0;
  }
  else { // We're checkpointing:
    outlog("Recovering from checkpoint %s\n", FLAGS_checkpoint_recover.c_str());
    
    if (read_checkpoint(FLAGS_checkpoint_recover.c_str(),
			model,     topic_suffstats, gammas,
			lhoods,    heldout_gammas,  heldout_lhoods,
			bound,     heldout_bound,   convergence,
			old_bound, iter,    	    final_iters_flag,
			(heldout != NULL),          last_iter, K) == 0) {
      outlog("Successfully recovered from checkpoint.\n** Starting from iteration %i **\n", iter);
      debug("CHECKR X2b: gammas[0,0]=%f\nOA",gsl_matrix_get(gammas,0,0));
      debug("CHECKR X2c: lhoods[0,0]=%f\n",gsl_matrix_get(lhoods,0,0));

	// AG: Maybe here we can dump doc sums
    }
    else {
      outlog("Failed!\nERROR: Checkpoint recovery failed, exiting...%s\n","");
      _exit(26);
    }
  }

  sprintf(root, "%s/lda-seq/", file_root);
  make_directory(root);

  sprintf(em_log_filename, "%s/em_log.dat", file_root);
  em_log = fopen(em_log_filename, "w");
  
  while (iter < FLAGS_lda_sequence_min_iter || ((final_iters_flag == 0 || 
	 convergence > LDA_SEQ_EM_THRESH) && iter <= FLAGS_lda_sequence_max_iter) && !last_iter) {
    
    if (!(iter < FLAGS_lda_sequence_min_iter || ((final_iters_flag == 0 || 
	  convergence > LDA_SEQ_EM_THRESH)  && iter <= FLAGS_lda_sequence_max_iter))) {
      last_iter = 1;
    }
    
    outlog("\nE iter %3d\n", iter);
    fprintf(em_log, "%17.14e %5.3e\n", bound, convergence);

    old_bound = bound;
    gsl_matrix_set_zero(gammas);
    gsl_matrix_set_zero(lhoods);

    if (heldout != NULL) {
      gsl_matrix_set_zero(heldout_gammas);
      gsl_matrix_set_zero(heldout_lhoods);
    }

    debug("DEBUG ES1%s\n","");
    for (k = 0; k < K; k++) 
      gsl_matrix_set_zero(topic_suffstats[k]);
    
    debug("DEBUG ES1a%s\n","");

    // compute the likelihood of a sequential corpus under an LDA seq model and find the evidence lower bound.
    bound = lda_seq_infer(model,  data, topic_suffstats, gammas, 
			  lhoods, iter, last_iter,       file_root, false); // Main part of the E-step
    
    if (heldout != NULL) {
      outlog("(Iteration %d): Fitting heldout SS-LM.\n",iter);
      heldout_bound = lda_seq_infer(model, heldout, NULL, heldout_gammas, 
				    heldout_lhoods, iter, last_iter, file_root, true);
    }

    debug("DEBUG ES2%s\n","");

    // print out the gammas and likelihoods.
    sprintf(name, "%s/gam-iter%d.dat", root, iter);
    mtx_fprintf(name, gammas);
    
    sprintf(name, "%s/lhoods-iter%d.dat", root, iter);
    mtx_fprintf(name, lhoods);
    
    if (heldout != NULL) {
      sprintf(name, "%s/heldout_lhoods-iter%d.dat", root, iter);
      mtx_fprintf(name, heldout_lhoods);

      sprintf(name, "%s/heldout_gam-iter%d.dat", root, iter);
      mtx_fprintf(name, heldout_gammas);
    }
    
    // fit the variational distribution
    // Main part of the M-step:
    double topic_bnd = fit_lda_seq_topics(model, topic_suffstats, iter); 
    bound           += topic_bnd;

    debug("DEBUG ES3%s\n","");    
    
    // If this is rDIM, update mu
    if (FLAGS_model == "rdim") {
      if (iter == 0) {
	outlog("%s", "Writing initial mu to file.\n");
	char name[400];
	sprintf(name, "%s/mu-init-K_by_nmetafields.dat", file_root);
	FILE* f = fopen(name, "w");
	gsl_matrix_fprintf(f, model->mu, "%.14f");
	fclose(f);
      }
      
      outlog("(Iteration %d): Updating mu...\n", iter);
      update_mu(data, model, iter); // Added by AG for rDIM
      
      outlog("%s", "done updating mu: writing to file.\n");
      char fname[400];
      sprintf(fname, "%s/mu-iter_%d-K_by_nmetafields.dat", file_root, iter);
      write_gsl_matrix_col_major(model->mu, fname);
    }
        
    write_lda_seq(model, root);

    if ((bound - old_bound) < 0) {
      if (LDA_INFERENCE_MAX_ITER == 1)
	LDA_INFERENCE_MAX_ITER = 2;
      if (LDA_INFERENCE_MAX_ITER == 2)
	LDA_INFERENCE_MAX_ITER = 5;
      if (LDA_INFERENCE_MAX_ITER == 5)
	LDA_INFERENCE_MAX_ITER = 10;
      if (LDA_INFERENCE_MAX_ITER == 10)
	LDA_INFERENCE_MAX_ITER = 20;

      outlog("**WARNING: Bound went down by %18.14f: Increasing var iterations to %d\n", 
	     bound-old_bound, LDA_INFERENCE_MAX_ITER);
    }

    // check for convergence
    convergence = fabs((bound - old_bound) / old_bound);

    if (convergence < LDA_SEQ_EM_THRESH) {
      final_iters_flag = 1;
      LDA_INFERENCE_MAX_ITER = 500;
      outlog("Starting final iterations: max iter = %d\n", LDA_INFERENCE_MAX_ITER);
      convergence = 1.0;
    }

    outlog("\n(%02d) LDA-SEQ bound=% 15.7f; convergence=% 15.7e\n", iter, bound, convergence);
    if (heldout != NULL) {
      outlog("Heldout bound=% 15.7f\n", heldout_bound);
    }

    iter++;

    // CHECKPOINT OUT:
    if (FLAGS_enable_checkpointing) {
      outlog("Exporting checkpoint for iter=%i\n",iter);
      if (write_checkpoint(FLAGS_checkpoint_outdir.c_str(), model,
			   topic_suffstats, gammas,lhoods, heldout_gammas, 
			   heldout_lhoods, bound, heldout_bound, convergence,
			   old_bound, iter, (heldout != NULL), 
			   final_iters_flag, last_iter, K) == 0) {
	outlog("Checkpoint successfully exported!\n%s","");
      }
      else {
	outlog("ERROR: Failed to export checkpoint at iteration %i, exiting...\n",iter);
	exit(-21);
      }
      if (FLAGS_die_after_checkpoint == true) {
        exit(21);
      }
    }
  }

  // added by YH
  for (k = 0; k < K; k++) 
    gsl_matrix_free(topic_suffstats[k]);
     
  gsl_matrix_free(gammas);
  gsl_matrix_free(lhoods);

  if (heldout != NULL) {
    gsl_matrix_free(heldout_gammas);
    gsl_matrix_free(heldout_lhoods);
  }
  
  return(bound);
}
 
#ifdef SPARSE
double g_round(double x, unsigned int digits) {
    double fac = pow(10, digits);
    return round(x*fac) / fac;
}

void compcol2dns(gsl_matrix* D, const gsl_spmatrix* C) {
  double v;

  for (size_t i = 0; i < C->size1; i++)
    for (size_t j = 0; j < C->size2; j++) {
      v = gsl_spmatrix_get(C,i,j);
      gsl_matrix_set(D,i,j,v);
    }
}

void add_noise_scaling(gsl_matrix *dest, const double epsilon) {
  assert(dest->size1 == dest->size2);

  double val;
  for (int i=0; i < dest->size1; i++) {
    val = gsl_matrix_get(dest,i,i);
    gsl_matrix_set(dest,i,i,val+epsilon);
  }
}

void gsl_spmatrix_get_row(gsl_vector* v, const gsl_spmatrix* M, const size_t x) {
  assert(x <= (M->size1 - 1) && v->size == M->size2);
  
  for (int i=0; i < M->size2; i++)
    gsl_vector_set(v, i, gsl_spmatrix_get(M,x,i));

  return;
}

void gsl_spmatrix_set_col(gsl_spmatrix* M, const size_t x, const gsl_vector* v) {
  for (int i=0; i < v->size; i++)
    gsl_spmatrix_set(M, i, x, gsl_vector_get(v, i));
}

void make_term1(const corpus_seq_t* data, lda_seq* model) {
  char hdf5fname[BUFSIZ];

  // We have an HDF version to read, and we're not recovering:
  if (strcmp("", FLAGS_read_hdf_term1.c_str()) != 0 && strcmp("", FLAGS_checkpoint_recover.c_str()) == 0) {
    model->eq32_term1 = read_hdf_term1(data->nmetafields);
    
    sprintf(hdf5fname, "%s/eq32_term1.csr3.h5", FLAGS_outname.c_str());
    
    hdf_write_gsl_to_csr3(model->eq32_term1,hdf5fname);
    outlog("Done writing A.%s\n","");
    gsl_matrix_free(model->eq32_term1);
  }
  else if (strcmp("", FLAGS_read_csr_term1.c_str()) != 0) {
    outlog("Expecting a CSR3 copy of eq_term1 to have been symlinked to the wd, not reading or creating.%s\n","");
  }
  else if (strcmp("", FLAGS_checkpoint_recover.c_str()) != 0) {
    // Hacky stuff here, but I want to avoid copying or reading the CSR3 matrix
    // It should be /somewhere/ near the checkpoint recovery location. (ie. two directories up)
    char fname[BUFSIZ];

    sprintf(hdf5fname, "%s/eq32_term1.csr3.h5", FLAGS_outname.c_str());
    sprintf(fname,     "%s/../../eq32_term1.csr3.h5", FLAGS_checkpoint_recover.c_str());

    outlog("Creating symlinks for recovery (hdf5fname='%s' fname='%s')\n",hdf5fname,fname);
    symlink(fname,hdf5fname);

    outlog("Expecting a CSR3 copy of eq_term1 in the checkpoint, not reading or creating.%s\n","");
  }
  else  {
    int nmetas = data->nmetafields;    

    gsl_spmatrix *t0 = gsl_spmatrix_alloc_nzmax(nmetas, nmetas,1,GSL_SPMATRIX_TRIPLET);
    gsl_spmatrix *S  = gsl_spmatrix_compcol(t0);
    gsl_spmatrix_free(t0);

    gsl_spmatrix_set_zero(S);

    t0 = gsl_spmatrix_alloc_nzmax(nmetas, nmetas,1,GSL_SPMATRIX_TRIPLET);
    gsl_spmatrix *eq32_term1_tmp2 = gsl_spmatrix_compcol(t0);
    gsl_spmatrix_free(t0);

#pragma omp parallel num_threads(FLAGS_threads)
    {
      gsl_spmatrix 
	*tmp             = gsl_spmatrix_alloc_nzmax(nmetas, nmetas,1,GSL_SPMATRIX_TRIPLET),
	*eq32_term1_tmp  = gsl_spmatrix_compcol(tmp),
	*eq32_term1_m    = gsl_spmatrix_alloc_nzmax(data->nmetafields,nmetas, 1, GSL_SPMATRIX_CCS),
	*tau_t_d         = gsl_spmatrix_alloc_nzmax(  nmetas, 1, 1, GSL_SPMATRIX_TRIPLET),
	*tau_t_d_trans   = gsl_spmatrix_alloc_nzmax(1,nmetas   , 1, GSL_SPMATRIX_TRIPLET),
	*runner_t;
      gsl_vector
	*tau_t_d_vec     = gsl_vector_alloc(nmetas);
      time_t start; 
      gsl_spmatrix_free(tmp);

      // We do this backwards because later time-steps tend to have more 
      // so doing them first will reduce straggling threads.
#pragma omp for schedule(dynamic)
      for (int t = data->len-1; t >= 0; t--) {
	tmp      = gsl_spmatrix_alloc_nzmax(nmetas, nmetas,1,GSL_SPMATRIX_TRIPLET);
	runner_t = gsl_spmatrix_compcol(tmp);
	gsl_spmatrix_free(tmp);
	
	gsl_spmatrix *t1, *t2;
	
	for (int d=0; d < data->corpus[t]->ndocs; ++d) {
	  gsl_vector_set_zero(tau_t_d_vec);
	  gsl_spmatrix_set_zero(tau_t_d);
	  gsl_spmatrix_set_zero(tau_t_d_trans);

	  gsl_spmatrix_get_row(tau_t_d_vec, data->tau[t], (size_t) d);
	  gsl_spmatrix_set_col(tau_t_d, (size_t) 0, tau_t_d_vec); // tau_t_d must be in TRIPLET format.
	  
	  gsl_spmatrix_transpose_memcpy(tau_t_d_trans, tau_t_d);
	  t1 = gsl_spmatrix_compcol(tau_t_d);
	  t2 = gsl_spmatrix_compcol(tau_t_d_trans);
	  
	  gsl_spblas_dgemm(1.0, t1, t2, eq32_term1_m);
	  
	  start=clock();
	  
	  gsl_spmatrix_memcpy(eq32_term1_tmp, runner_t);
	  gsl_spmatrix_add(runner_t, eq32_term1_tmp, eq32_term1_m);
	  
	  outlog("Building Eq 32 term1 (t=%d / %d, d=%d / %d ) took time = %f\n", 
		 t+1, (int) data->len, d+1, (int) data->corpus[t]->ndocs, (double) (clock() - start));
	  
	  gsl_spmatrix_free(t1);
	  gsl_spmatrix_free(t2);
	}
	
#pragma omp critical
	{
	  gsl_spmatrix_memcpy(eq32_term1_tmp2, S);
	  gsl_spmatrix_add(S, eq32_term1_tmp2, runner_t);
	}
	
	gsl_spmatrix_free(runner_t);
      }
      
      gsl_spmatrix_free(tau_t_d);
      gsl_spmatrix_free(tau_t_d_trans);
      gsl_vector_free(tau_t_d_vec);
      gsl_spmatrix_free(eq32_term1_tmp);
      gsl_spmatrix_free(eq32_term1_m);
    } // End GOMP pragma
    
    model->eq32_term1 = gsl_matrix_calloc(nmetas,nmetas);

    outlog("Incorporating hyper-prior for \\mu into term1 (ie. estimating as \\hat{\\mu}\n%s","");
    double term1_diag_val = 1.0;

    if (FLAGS_sigma_l != FLAGS_sigma_mu)
      term1_diag_val = (FLAGS_sigma_l * FLAGS_sigma_l) / (FLAGS_sigma_mu * FLAGS_sigma_mu);

#pragma omp parallel shared(S, model)
    {
#pragma omp for schedule(static)
      for (size_t i = 0; i < S->size1; i++) {
	if ((int) i % 1000 == 0) {
	  outlog("Densifying term 1 (row %i of %i)\n", (int) i, (int) S->size1);
	}
	gsl_matrix_set(model->eq32_term1, i, i, gsl_spmatrix_get(S, i, i) + term1_diag_val);
	
	for (size_t j = i+1; j < S->size2; j++)
	  gsl_matrix_set(model->eq32_term1, i, j,gsl_spmatrix_get(S, i, j));
      }
    }
    
    int nz;
#pragma omp parallel shared(nz)
    {
#pragma omp for reduction(+:nz) schedule(static)
      for (long i=0; i < nmetas*nmetas; i++)
	if (model->eq32_term1->data[i] != 0.0)
	  nz++;
    }//GOMP
    model->eq32_term1_nz = nz;
    
    // Broke down and now we require HDF5:
    outlog("Writing eq_term1 to the HDF5 CSR3 file '%s'...\n",FLAGS_outname.c_str());
    sprintf(hdf5fname, "%s/eq32_term1.csr3.h5", FLAGS_outname.c_str());

    hdf_write_gsl_to_csr3(model->eq32_term1,hdf5fname);

    //outlog("Writing eq_term1 to the HDF5 file '%s'...\n",FLAGS_outname.c_str());
     
    //sprintf(hdf5fname, "%s/eq32_term1.h5", FLAGS_outname.c_str());
    //hdf_write_matrix(model->eq32_term1, hdf5fname);

    outlog("Done writing A.%s\n","");
    gsl_matrix_free(model->eq32_term1);
  }

  return;
}

// Updating mu (the rDIM metadata coefficients)
// Equation 32 in YH's derivation
// Added by AG
void update_mu(const corpus_seq_t* data, lda_seq* model, const int iter) { // SPARSE
  gsl_matrix *rhs = gsl_matrix_alloc(data->nmetafields, model->ntopics);

#pragma omp parallel num_threads(FLAGS_kthreads)
{
  double l_t_d_k, ival; // Influence of document d, and time t, in topic k.
  gsl_vector
    *eq32_term2  = gsl_vector_alloc(data->nmetafields),
    *tau_t_d_vec = gsl_vector_alloc(data->nmetafields);

#pragma omp for
  for (int k=0; k < model->ntopics; ++k) { // For every topic, k
    gsl_vector_set_zero(eq32_term2);
    
    outlog("Updating mu (k=%i)\n",k);
    for (int t=0; t < data->len; ++t) {
      for (int d=0; d < data->corpus[t]->ndocs; ++d) {
	gsl_spmatrix_get_row(tau_t_d_vec, data->tau[t], (size_t) d);
	l_t_d_k = g_round(gsl_matrix_get(model->influence->renormalized_doc_weights[t],d,k), PREC);
	gsl_vector_scale(tau_t_d_vec, l_t_d_k);
	gsl_vector_add(eq32_term2, tau_t_d_vec);
      }
    }
    
#ifdef MKL
#pragma omp critical(rhs_setter) 
    {
      gsl_matrix_set_col(rhs, k, eq32_term2);
    }
#else
    // This non-MKL sparse solver is garbage // use MKL!
    if (SparseSolve(model->eq32_term1, eq32_term2, &new_mu_k.vector, 3 * sqrt(model->eq32_term1->size1)) != GSL_SUCCESS) {
      outlog("FATAL FAILURE: Failed to update mu because term 1 is degenerate: exiting.%s\n","");
      exit(-2);
    }
    else {
      outlog("SUCCESS: Solver converged! Hurray!%s):\n","");
    }
#endif
  }

  gsl_vector_free(eq32_term2);
  gsl_vector_free(tau_t_d_vec);

  } // End pragma

#ifdef MKL
  outlog("Calling MKL::LAPACKE_dgesvx() for multple right-hand sides linear solve...\n%s","");
 
  // MPI NOTE: This function is a large, dense linear solve from MKL
  // There are MKL calls that use MPI for this (I think PARDISO and ScaLAPACK)
  // If MPI elsewhere is too hard, this would be great to have MPI'd

  char A_fname[BUFSIZ]; 
  sprintf(A_fname, "%s/eq32_term1.csr3.h5", FLAGS_outname.c_str()); 
  int result = fork_pardiso_solve(A_fname, rhs, model->mu, model->eq32_term1_nz); 
  //int result = MKLSolve_lapacke(model->eq32_term1, rhs, model->mu, model->A_factor, model->try_factor, model->A_factored); 
  //int  result = MKLSolve_pardiso(model->eq32_term1, rhs, model->mu, model->eq32_term1_nz);

  if (result != 0) { 
    outlog("FAILURE: PARDSIO failed with error=%i (see stdout for PARDISO report)\n",result); 
    _exit(-2); 

    outlog("Trying MKL P-Direct Solver...%s\n",""); 

    result = MKLSolve_lapacke(model->eq32_term1, rhs, model->mu, model->A_factor, model->try_factor, model->A_factored); 
    if (result != 0) {     
      outlog("FATAL FAILURE: Failed to update mu (see output from mkl_solve): exiting.%s\n",""); 
      _exit(-2); 
    } 
  }  

  gsl_matrix_free(rhs);

#endif
  
  debug("mu[0,0]=%f\nmu[0,1]=%f\nmu[1,1]=%f\nmu[-1,-1]=%f\n",
	gsl_matrix_get(model->mu, 0, 0),
	gsl_matrix_get(model->mu, 0, 1),
	gsl_matrix_get(model->mu, 1, 1),
	gsl_matrix_get(model->mu, data->nmetafields-1, model->ntopics-1));
}

#endif

/*
 * read and write lda sequence variational distribution
 *
 */
void write_lda_seq(const lda_seq* model, const char* root) {
  char name[400];
  sprintf(name, "%sinfo.dat", root);
  FILE* f = fopen(name, "w");

  params_write_int(f, (char* const)   "NUM_TOPICS", model->ntopics);
  params_write_int(f, (char* const)   "NUM_TERMS",  model->nterms);
  params_write_int(f, (char* const)   "SEQ_LENGTH", model->nseq);
  params_write_gsl_vector(f, (char* const) "ALPHA", model->alpha);
  
  fclose(f);

  for (int k = 0; k < model->ntopics; k++) {
    const int tmp = k;
    outlog("\nwriting topic %03d\n", tmp);
    sprintf(name, "%stopic-%03d", root, tmp);
    write_sslm_var(model->topic[tmp], name);
  }
  
  if (FLAGS_model == "fixed" || FLAGS_model == "rdim") {
    for (int t=0; t < model->influence->ntime; ++t) {
      sprintf(name, "%sinfluence_time-%03d", root, t);
      outlog("\nwriting influence weights for time %d to %s\n", t, name);
      gsl_matrix* influence_t = model->influence->doc_weights[t];
      assert(model->ntopics == influence_t->size2);
      mtx_fprintf(name, influence_t);

      sprintf(name, "%srenormalized_influence_time-%03d", root, t);
      outlog("\nwriting influence weights for time %d to %s\n", t, name);
      influence_t = model->influence->renormalized_doc_weights[t];
      assert(model->ntopics == influence_t->size2);
      mtx_fprintf(name, influence_t);
    }
  }
}

// Read information about a particular model.
// This model should be named "{root}info.dat"
// and should contain the following rows:
// number_topics
// number_times
// alpha, as a gsl vector
lda_seq* read_lda_seq(const char* root, corpus_seq_t* data) {
  char name[400];
  lda_seq* model = (lda_seq*) malloc(sizeof(lda_seq));
  
  sprintf(name, "%sinfo.dat", root);
  FILE* f = fopen(name, "r");
  
  if (f == NULL) {
    outlog("Unable to open file %s.  Failing.\n", name);
    exit(1);
  }

  params_read_int(f, (char* const)   "NUM_TOPICS", &(model->ntopics));
  params_read_int(f, (char* const)    "NUM_TERMS", &(model->nterms));
  params_read_int(f, (char* const)   "SEQ_LENGTH", &(model->nseq));
  params_read_gsl_vector(f, (char* const) "ALPHA", &(model->alpha));
  fclose(f);

  model->topic = (sslm_var**)  malloc(sizeof(sslm_var*)   * model->ntopics);
  // TODO (long-term) AG: write out mus in write_lda_seq 
  // TODO (long-term) AG: read in mus here.

  for (int k = 0; k < model->ntopics; k++) {
    outlog( "Reading topic %d\n", k);
    sprintf(name, "%stopic-%03d", root, k);
    model->topic[k] = read_sslm_var(name);
    
    model->topic[k]->w_phi_l    = gsl_matrix_alloc(model->nterms, model->nseq);
    model->topic[k]->w_phi_sum  = gsl_matrix_alloc(model->nterms, model->nseq);
    model->topic[k]->w_phi_l_sq = gsl_matrix_alloc(model->nterms, model->nseq);
    
    if (FLAGS_model == "dim") {
      sprintf(name, "%sw_phi_l-%d", root, k);
      mtx_fscanf(name, model->topic[k]->w_phi_l);
      
      sprintf(name, "%sw_phi_sum-%d", root, k);
      mtx_fscanf(name, model->topic[k]->w_phi_sum);
      
      sprintf(name, "%sw_phi_l_sq-%d", root, k);
      mtx_fscanf(name, model->topic[k]->w_phi_l_sq);
    }
  }
  
  if (FLAGS_model == "dim" || FLAGS_model == "rdim" && data != NULL) {
    model->influence = (inf_var*) malloc(sizeof(inf_var));
    model->influence->doc_weights = (gsl_matrix**) malloc(sizeof(gsl_matrix*));    
    int t;
    model->influence->ntime = model->nseq;

    for (t=0; t < model->nseq; ++t) {
      sprintf(name, "%sinfluence_time-%03d", root, t);
      outlog("\nReading influence weights for time %d from %s\n", t, name);
      model->influence->doc_weights[t] = gsl_matrix_alloc(data->corpus[t]->ndocs, model->ntopics);
      mtx_fscanf(name, model->influence->doc_weights[t]);
    }
  } else
    model->influence = NULL;
  
  return(model);
}

/*
 * update lda sequence sufficient statistics from an lda posterior
 *
 */
void update_lda_seq_ss(int time, const doc_t* doc, const lda_post* post, gsl_matrix** ss) {
  int K = post->phi->size2, N = doc->nterms;

#pragma omp parallel num_threads(FLAGS_kthreads)
{
  gsl_matrix* topic_ss;

#pragma omp for schedule(static)
  for (int k = 0; k < K; k++) {
    topic_ss = ss[k];
    for (int n = 0; n < N; n++)
      minc(topic_ss, doc->word[n], time, doc->count[n] * mget(post->phi, n, k));
    }
  }// GOMP
}


void update_lda_seq_ss_serial(int time, const doc_t* doc, const lda_post* post, gsl_matrix** ss) {
  int K = post->phi->size2, N = doc->nterms;

  gsl_matrix* topic_ss;

  for (int k = 0; k < K; k++) {
    topic_ss = ss[k];
    for (int n = 0; n < N; n++)
      minc(topic_ss, doc->word[n], time, doc->count[n] * mget(post->phi, n, k));
  }

}

/*
 * fit lda sequence
 *
 */
double fit_lda_seq_topics(lda_seq* model, gsl_matrix** ss, const int iter) {
  double lhood = 0;
  outlog("(Iteration %d) M Step.\n", iter);
  //set_thread_affinity();
  //omp_set_max_active_levels(2);
  //omp_set_nested(1);
  //mklspace::mkl_set_dynamic(0);
  //omp_set_dynamic(0);

#pragma omp parallel num_threads(FLAGS_kthreads)
  {
#pragma omp for reduction(+:lhood) schedule(static)
    for (int k = 0; k < model->ntopics; k++) { // Line 13 of rDIM algorithm
      outlog("(Iteration %d): fitting topic %d (thread %d)\n", 
	     iter, k, (int) omp_get_thread_num());
      double val = fit_sslm(model->topic[k], model, ss[k]); 
      lhood     += val;
    }
  } // End Pragma
  
  return(lhood);
}

/*
 * allocate lda seq
 *
 */
lda_seq* new_lda_seq(corpus_seq_t* data, int W, int T, int K) {
  lda_seq* model = (lda_seq*) malloc(sizeof(lda_seq));
  
  model->ntopics = K;
  model->nterms  = W;
  model->nseq    = T;
  model->alpha   = gsl_vector_alloc(K);
  model->topic   = (sslm_var**)  malloc(sizeof(sslm_var*) * K);
  
  model->eq32_term1_nz = 0;
  model->A_factored    = false; // Have not solved yet -- after first solve for mu, this is true
  model->try_factor    = true;
  model->A_factor      = (mkl_factor*) malloc(sizeof(mkl_factor));
  model->A_factor->n   = 0;

  if (FLAGS_model == "rdim") {
    model->mu = init_mu(data->nmetafields, K); // Now COL_MAJOR
  }

  // Create the vectors of total counts for each time.
  model->influence_sum_lgl = (gsl_matrix**) malloc(sizeof(gsl_matrix*) * K);
  
  for (int k = 0; k < K; k++) {
    model->influence_sum_lgl[k] = gsl_matrix_calloc(W, T);
    model->topic[k] = sslm_var_alloc(W, T);

    if (k < FLAGS_fix_topics) 
      model->topic[k]->chain_variance = 1e-10;
      
    model->topic[k]->w_phi_sum  = gsl_matrix_calloc(W, T);
    model->topic[k]->w_phi_l_sq = gsl_matrix_calloc(W, T);
  }
  
  model->influence = inf_var_alloc(K, data);

  return(model);
}

// added by YH
void free_lda_seq_model(lda_seq* m) {
  int K = m->ntopics;

  for (int k = 0; k < K; k++) {
    gsl_matrix_free(m->influence_sum_lgl[k]);
    // gsl_matrix_free(m->topic[k]->w_phi_l);
    gsl_matrix_free(m->topic[k]->w_phi_sum);
    gsl_matrix_free(m->topic[k]->w_phi_l_sq);
    sslm_var_free(m->topic[k]);
  }

  free(m->influence_sum_lgl);
  free(m->topic);
  gsl_vector_free(m->alpha);
  inf_var_free(m->influence);

  if (FLAGS_model == "rdim") {
    gsl_matrix_free(m->mu); // Added by AG
    //    gsl_matrix_free(m->eq32_term1);
  }

  free(m);
}

/*
 * initialize from sufficient statistics (expected counts).
 *
 */
void init_lda_seq_from_ss(lda_seq* model, double topic_chain_variance, double topic_obs_variance, double alpha, gsl_matrix* init_suffstats) {
    gsl_vector_set_all(model->alpha, alpha);

    for (int k = 0; k < model->ntopics; k++) {
        gsl_vector slice = gsl_matrix_column(init_suffstats, k).vector;
        sslm_counts_init(model->topic[k], topic_obs_variance, topic_chain_variance, &slice);

	if (k < FLAGS_fix_topics) 
	  model->topic[k]->chain_variance = 1e-10;
	
	// commented by YH
	//model->topic[k]->w_phi_l    = gsl_matrix_calloc(model->nterms, model->nseq);
	//model->topic[k]->w_phi_sum  = gsl_matrix_calloc(model->nterms, model->nseq);
	//model->topic[k]->w_phi_l_sq = gsl_matrix_calloc(model->nterms, model->nseq);
    }
}

// Initializ mu (the coefficients in our meta-data regression)
// Currently this is done by sampling form multivariate Gaussians for every topic
// The assumption here (currently) is that the metadata are independently random
// with respect to topic.
// The values, drawn from N(0,FLAGS_sigma_mu), or normalized to [0,1] after being drawn.
//
// Added by AG
gsl_matrix* init_mu(const int nmetafields, const int ntopics) {
  gsl_matrix *M      = gsl_matrix_alloc(nmetafields, ntopics); // Number of docs in corpus[t]
  gsl_vector* values = gsl_vector_alloc(ntopics);
  double      max, min, val, sigma = FLAGS_sigma_mu;
  
  outlog("Initializing mu (%d x %d) from %d multivariate Gaussians of the form: N(0,%g).\n", 
	 ntopics, nmetafields, ntopics, sigma);

  // For each column (ie. field / variable) of metadata:
  for (int s=0; s < nmetafields; ++s) {
    gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);

    gsl_vector_set_zero(values);

    // This is important to maintain determinism when setting --rng_seed
    // while retaining randomness across multivariates:
    gsl_rng_set(r, (long) FLAGS_rng_seed + (long) s); // Must be addition...

    // Drawing the random values:
    for (int k=0; k < ntopics; ++k) {
      val    = gsl_ran_gaussian(r, sigma);
      gsl_vector_set(values, (size_t) k, val);
    }
    
    // If k == 1, then min == max, and we shouldn't norm:
    // If k == 2, then max-min == 1.0 and all values will be in {0,1}...would rather have unnormed doubles.
    if (ntopics > 2) { 
      // For norming:
      min = gsl_vector_min(values);
      max = gsl_vector_max(values);
      assert(max != min);

      // Norming:
      for (int k=0; k < ntopics; ++k) {
	val = gsl_vector_get(values, (size_t) k);
	val = (val-min) / (max-min);
	//	debug("DEBUG: mu_norm2 (val=%f)\n", val);
	gsl_vector_set(values, (size_t) k, val);
      }
    }

    // Setting:
    gsl_matrix_set_row(M, (size_t) s, values); // Used to set col, now row -- for COL_MAJOR linear solve
    gsl_rng_free(r);
  }
  
  gsl_vector_free(values);

  return M;
}

// added by YH
void free_corpus(corpus_t* data) {
  if (data == NULL)
    return;

  for (int i = 0; i < data->ndocs; i++) {
    free(data->doc[i]->word);
    free(data->doc[i]->count);
    free(data->doc[i]->lambda);
    free(data->doc[i]->log_likelihoods);
    free(data->doc[i]);
  }

  free(data->doc);
  free(data);
  data = NULL;
}

// added by YH
void free_corpus_seq(corpus_seq_t* data_seq) {
  if (data_seq == NULL)
    return;

  outlog("%s","Freeing corpus...");

  for (int i = 0; i < data_seq->len; i++) {
    for (int j = 0; j < data_seq->corpus[i]->ndocs; j++) {
      free(data_seq->corpus[i]->doc[j]->word);
      free(data_seq->corpus[i]->doc[j]->count);
      free(data_seq->corpus[i]->doc[j]->lambda);
      free(data_seq->corpus[i]->doc[j]->log_likelihoods);
      free(data_seq->corpus[i]->doc[j]);
    }

    // Added by AG
    if (FLAGS_model == "rdim") {
#ifdef SPARSE
      gsl_spmatrix_free(data_seq->tau[i]);
#else
      gsl_matrix_free(data_seq->tau[i]);
#endif
    }
    free(data_seq->corpus[i]->doc);
    free(data_seq->corpus[i]);
  }
  // Added by AG
  if (FLAGS_model == "rdim")
    free(data_seq->tau);

  free(data_seq->corpus);
  free(data_seq);

  data_seq = NULL;
  outlog("%s","done.\n");
}

// Note: Does not free(ptr)
void mkl_spmat_free(mkl_spmat* ptr) {
  mklspace::mkl_free(ptr->Mv);
  mklspace::mkl_free(ptr->Mi);
  mklspace::mkl_free(ptr->Mj);
}
