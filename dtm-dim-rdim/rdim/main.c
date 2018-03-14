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

#include <stdio.h>
#include <error.h>
#include <unistd.h>
#include <signal.h>
#include <fenv.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sched.h>

#include "gflags.h"
#include "main.h"
#include "gsl/gsl_ieee_utils.h"
#include "hdf5.h"
#include "hdf5_hl.h"

namespace mklspace
{
#include "mkl_lapacke.h"
#include "mkl.h"
#include "mkl_types.h"
#include "mkl_spblas.h"
#include "mkl_service.h"
}

#ifdef SPARSE
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>
#include <gsl/gsl_spblas.h>
#include <omp.h>
#endif

#include "lda-seq.h"

#ifdef CUDA
#include <cuda_runtime.h>
#include <cublas.h>
#include <cusparse.h>

cusparseHandle_t hndl;
#endif

extern int parallelism_enabled;

#define PREC 15 
#define PARENT_WRITE_PIPE  0
#define PARENT_READ_PIPE   1
#define READ_FD  0
#define WRITE_FD 1
#define PARENT_READ_FD  ( pipes[PARENT_READ_PIPE][READ_FD]   )
#define PARENT_WRITE_FD ( pipes[PARENT_WRITE_PIPE][WRITE_FD] ) 
#define CHILD_READ_FD   ( pipes[PARENT_WRITE_PIPE][READ_FD]  )
#define CHILD_WRITE_FD  ( pipes[PARENT_READ_PIPE][WRITE_FD]  )

DEFINE_string(mode,             "fit", "The function to perform. " "Can be fit, est, or time.");
DEFINE_string(model,            "dtm", "The function to perform. " "Can be dtm or dim.");
DEFINE_string(corpus_prefix,    "",    "The prefix of the corpus to use.");
DEFINE_string(lda_model_prefix, "",    "The name of a fit model to be " "used for testing likelihood.  Appending \"info.dat\" " "to this should give the name of the file.");
DEFINE_int32(heldout_time,     -1,     "A time up to (but not including) which we wish to train, " "and at which we wish to test.");
DEFINE_string(output_table,     "",    "");
DEFINE_string(params_file,      "settings.txt", "A file containing parameters for this run.");
DEFINE_bool(initialize_lda,     false,          "If true, initialize the model with lda.");

DEFINE_string(outname,          "",    "");
DEFINE_double(top_obs_var,      0.5,   "");
DEFINE_double(top_chain_var,    0.005, "");
DEFINE_double(alpha,           -10.0,  "");
DEFINE_double(ntopics,         -1.0,   "");
DEFINE_int32(lda_max_em_iter,   20,    "");
DEFINE_string(heldout_corpus_prefix,   "", "");
DEFINE_string(initial_lda_ss,   "",    "Path to a sufficient stats matrix for the initializing LDA");
DEFINE_int32(metafields,        0,     "How many fields are apparent in the metadata files (\tau_t). All fields must be present in each time-slice.");
DEFINE_int32(debug,             0,    "Print debug output? Not for the faint of heart...");
DEFINE_int32(kthreads,          1,    "The number of threads for K-indexed loops");
DEFINE_int32(threads,           1,    "The number of threads");

extern int LDA_INFERENCE_MAX_ITER;

// If OpenBLAS was compiled without GNU-compliance enabled, threads can jump cores...
void set_thread_affinity() {
  if (FLAGS_threads > 1) {
#pragma omp parallel num_threads(FLAGS_threads)
  { 
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(omp_get_thread_num(), &mask);
    int res = sched_setaffinity(0, sizeof(mask), &mask);
    
    if (res == -1) {
      outlog("Failed set_thread_affinity(): Non-GNU system?\n%s","");
      outlog("Would normally bail, but trying anyway (might have competing GOMP threads%s\n","");
      //_exit(-6);
    }
  }
}
}

// SSE 4.5 enabled CPUs:
void set_fpu (unsigned int mode) {
  asm("fldcw %0" : : "m" (*&mode));
}

/*
 * fit a model from data
 *
 */
void fit_dtm(int min_time, int max_time) {
    char name[400];
    char run_dir[400];
    sprintf(run_dir, "%s/", FLAGS_outname.c_str());
    if (!directory_exist(run_dir))
      make_directory(run_dir);

    gsl_matrix* topics_ss;
    corpus_t* initial_lda_data = read_corpus(FLAGS_corpus_prefix.c_str());
    if (FLAGS_initialize_lda) {
      // initialize (a few iterations of LDA fitting)
      outlog("%s","Initializing dynamic model from a static LDA\n");
      outlog("Data: %s\n", FLAGS_corpus_prefix.c_str());
      
      lda* lda_model = new_lda_model(FLAGS_ntopics, initial_lda_data->nterms);
      gsl_vector_set_all(lda_model->alpha, FLAGS_alpha);
      
      lda_suff_stats* lda_ss = new_lda_suff_stats(lda_model);
      
      // initialize_lda_ss_from_data(initial_lda_data, lda_ss);
      initialize_lda_ss_from_random(initial_lda_data, lda_ss);
      
      // sgerrish: Why do we only define the topics once?
      lda_m_step(lda_model, lda_ss);
      
      sprintf(name, "%s/initial-lda", run_dir);

      LDA_INFERENCE_MAX_ITER = 25;

      lda_em(lda_model, lda_ss, initial_lda_data, FLAGS_lda_max_em_iter, name);
      sprintf(name, "%s/initial-lda-ss.dat", run_dir);
       
      write_lda_suff_stats(lda_ss, name);
      topics_ss = lda_ss->topics_ss;

      free(lda_ss);              // added by YH
      lda_ss = NULL;             // added by YH
      free_lda_model(lda_model); // added by YH
    } 
    else if (strcmp("",FLAGS_initial_lda_ss.c_str())) {
      outlog("Loading initial LDA  from %s (%i terms)...\n", 
	     FLAGS_initial_lda_ss.c_str(), initial_lda_data->nterms);
      
      topics_ss = gsl_matrix_calloc(initial_lda_data->nterms, FLAGS_ntopics);
      sprintf(name, "%s", FLAGS_initial_lda_ss.c_str());
      
      mtx_fscanf(name, topics_ss);
    }
    else {
      outlog("Need to specify either --initiali_lda=true or --initial_lda_ss\n%s","");
      outlog("No initializing LDA specified (either to fit, or by a file): exiting...\n%s","");
      _exit(-12);
    }

    free_corpus(initial_lda_data); // added by YH

    // estimate dynamic topic model
    corpus_seq_t
      *data_subset,
      *data_full    = read_corpus_seq(FLAGS_corpus_prefix.c_str(), FLAGS_ntopics, FLAGS_metafields),
      *heldout_data = NULL; // Added by AG

    // If we got a hold-out (ie. test) corpus:
    if (strcmp("",FLAGS_heldout_corpus_prefix.c_str()) != 0) {
      outlog("Holdout (ie. test) corpus specified: '%s'\n", FLAGS_heldout_corpus_prefix.c_str());
      heldout_data = read_corpus_seq(FLAGS_heldout_corpus_prefix.c_str(), FLAGS_ntopics, FLAGS_metafields);
    }

    if (max_time >= 0) { // We are training on a subset of the data.
      outlog("Caught max_time=%d and min_time=%d, preparing data_subset...", max_time, min_time);
      assert(max_time > min_time && min_time >= 0 && max_time < data_full->len);

      data_subset         = (corpus_seq_t*) malloc(sizeof(corpus_seq_t));
      //      data_subset->len    =  max_time - min_time + 1; // This may be wrong // AG
      data_subset->len    =  max_time - min_time; 
      data_subset->nterms =  data_full->nterms;
      data_subset->corpus = (corpus_t**) malloc(sizeof(corpus_t*) * data_subset->len);

      int max_nterms = 0, 
	       ndocs = 0;
       
      for (int i = min_time; i < max_time; ++i) {
	corpus_t* corpus = data_full->corpus[i];

	if (max_nterms > corpus->nterms)
	  max_nterms = max_nterms;
	else 
	  max_nterms = corpus->nterms;

	data_subset->corpus[i - min_time] = corpus;
	ndocs += corpus->ndocs;

	if (FLAGS_model == "rdim") {
#ifdef SPARSE
	  gsl_spmatrix* tau = data_full->tau[i];
#else
	  gsl_matrix*   tau = data_full->tau[i];
#endif
	  data_subset->tau[i-min_time] = tau;
	}
      }

      data_subset->max_nterms  = max_nterms;
      data_subset->ndocs       = ndocs;

      if (FLAGS_model == "rdim")
	data_subset->nmetafields = data_full->nmetafields; // Added by AG

      outlog("%s\n","done.");
    } 
    else { // Use the entire dataset.
      data_subset = data_full;
    }

    debug("main.c: Creating an LDA init...\n%s","");
    // Fit an LDA model to the whole (not sliced) data
    lda_seq* model_seq = new_lda_seq(data_subset, data_subset->nterms, data_subset->len, FLAGS_ntopics);

    debug("(L-6) model_seq->mu[k=0][0] = %.15f\n", gsl_matrix_get(model_seq->mu, 0, 0));

    // Initialize a DTM with the above LDA model:
    debug("main.c: Calling init_lda_seq_from_ss():\n%s","");

    init_lda_seq_from_ss(model_seq, FLAGS_top_chain_var, FLAGS_top_obs_var, FLAGS_alpha, topics_ss);

    // Run DTM / DIM / rDIM:
    debug("Fitting the SS-LM:\n%s","");
    
    fit_lda_seq(model_seq, data_subset, heldout_data, run_dir);

    gsl_matrix_free(topics_ss);      // added by YH

    if (max_time < 0) {
      free_corpus_seq(data_full);    // added by YH
      data_subset = NULL;            // added by YH
      free_lda_seq_model(model_seq); // added by YH

      return;
    }

    outlog("Computing document held-out LL at time %d", max_time);
    
    // Now find the posterior likelihood of the next time slice
    // using the most-recently-known time slice.
    lda* lda_model = new_lda_model(model_seq->ntopics, model_seq->nterms);
    make_lda_from_seq_slice(lda_model, model_seq, max_time - 1);

    int max_nterms = compute_max_nterms(data_full);

    // added by YH
    lda_post* post   = new_lda_post(model_seq->ntopics, max_nterms);
    post->model      = lda_model;
    post->doc_weight = NULL;
    double* table    = (double*) malloc(sizeof(double) * data_full->corpus[max_time]->ndocs);

    for (int d = 0; d < data_full->corpus[max_time]->ndocs; d++) {
      // added by YH
      post->doc = data_full->corpus[max_time]->doc[d];
      table[d]  = fit_lda_post(d, max_time, post, NULL, NULL, NULL, NULL, NULL);
    }
    
    char tmp_string[400];
    sprintf(tmp_string, "%s/heldout_post_%d.dat", FLAGS_outname.c_str(), max_time);
    FILE* post_file = fopen(tmp_string, "w");
    
    for (int d = 0; d < data_full->corpus[max_time]->ndocs; ++d)
      fprintf(post_file, "%f\n", table[d]);
    
    // Added by YH:
    fclose(post_file);
    free_corpus_seq(data_full);
    data_subset = NULL; 
    free(table); 
    free_lda_model(lda_model);
    free_lda_post(post);
}

#ifdef GTESTMODE
void mkl_spblas_dgemm_main(gsl_matrix *A, const gsl_matrix* B, gsl_matrix* C) {
  assert(C->size1 == C->size2);
  assert(A->size2 == C->size2);
  assert(B->size2 == C->size2);
  assert(A->size1*A->size2 == B->size1*B->size2);

  debug("SPMKL 1%s\n",""); 

  MKL_INT job[6]   = {0, 0, 1, 2, NULL, 1}; 
  // {dense->scr, input is 0-indexed, output should be 1-indexed, 
  //  translate entire matrix, N_nonzeros in each matrix (set after counting), 
  //  generate Xscr, xi and xj as output}

  const MKL_INT an = A->size1;
  const MKL_INT am = A->size2;
  const MKL_INT bn = B->size1;
  const MKL_INT bm = B->size2;
  MKL_INT      ldc = C->size1; // /MUST/ be copied.
  MKL_INT       NZ = 0;
  MKL_INT* AI = (MKL_INT*) malloc(sizeof(MKL_INT) * (an+1));
  MKL_INT* BI = (MKL_INT*) malloc(sizeof(MKL_INT) * (bn+1));
  MKL_INT     info;       // Output status from MKL
  const char trans = 'T'; // Transpose the 'A' argument to MKL. 
                    // Here we transpose B because the col/row-major translation transpose A and C.
  debug("SPMKL 2%s\n",""); 
  for (int i = 0; i < an*am; i++)
    if (A->data[i] != 0.0) NZ++;
  debug("SPMKL 3 NZ=%i\n",NZ); 

  double* Acsr = (double*) malloc(sizeof(double) * NZ);
  double* Bcsr = (double*) malloc(sizeof(double) * NZ);
  MKL_INT*  AJ = (MKL_INT*) malloc(sizeof(MKL_INT) * NZ);
  MKL_INT*  BJ = (MKL_INT*) malloc(sizeof(MKL_INT) * NZ);

  job[4] = NZ ;
  mklspace::mkl_ddnscsr(job, &an, &am, A->data, &am, Acsr, AJ, AI, &info);
  mklspace::mkl_ddnscsr(job, &bn, &bm, B->data, &bm, Bcsr, BJ, BI, &info);

  omp_set_num_threads(1);
  mklspace::mkl_set_num_threads(FLAGS_threads);
  mklspace::mkl_dcsrmultd(&trans, &bn, &bm, &am, Bcsr, BJ, BI, Acsr, AJ, AI, C->data, &ldc);
  mklspace::mkl_set_num_threads(1);
  omp_set_num_threads(FLAGS_threads); 

  debug("SPMKL 7%s\n","");
  free(Acsr);
  free(Bcsr);
  free(AJ);
  free(BJ);
  free(AI);
  free(BI);
  debug("SPMKL 8%s\n","");

  return;
}

int MKLSolve_lapacke_single_main(gsl_matrix* gsl_A, gsl_vector* gsl_b, gsl_vector* gsl_x) {
  MKL_INT
    n    = gsl_A->size1,
    nrhs = 1,
    info, ipiv[n]; // ?Or nrhs?  
  double  rcond, berr[nrhs], ferr[nrhs], rpivot;
  char equed = 'N';

  info = mklspace::LAPACKE_dgesvx(LAPACK_ROW_MAJOR,   'N', 'N',    n,  nrhs, gsl_A->data,  
				  n,    NULL,  n, ipiv, &equed, NULL, NULL,  gsl_b->data,  
				  nrhs, gsl_x->data, nrhs, &rcond, ferr, berr,   &rpivot);        
  if (info == (n+1)) {
    fprintf(stderr, "(LAPACKE linear solve) WARNING: Condition number below MACHINE_EPSILON.\n", (int) info);
  }
  else if (info > 0) {
    fprintf(stderr, "(LAPACKE linear solve) Matrix looks degenerate (retcode %i).\n", (int) info);
  }
  else if (info < 0) {
    outlog("(LAPACKE linear solve) threw an error: %i, exiting...\n", (int) info);
    _exit(info);
  }

  return(info);
}

void hdf_write_matrix_main(const gsl_matrix* M, const char* fname) {
  hid_t fid;
  hsize_t dims[2] = {M->size1, M->size2};
  herr_t status;

  printf("In hdf_write_matrix_main(): writing...\n");
  fid    = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  status = H5LTmake_dataset(fid, "term1", 2, dims, H5T_NATIVE_DOUBLE, M->data);
  assert(status != -1);

  status = H5Fclose(fid); assert(status != -1);
  assert(status != -1);

  return;
}

void hdf_write_two_matrix_main(const gsl_matrix* M1,const gsl_matrix* M2, const int att1, const char* fname) {
  hid_t fid;
  hsize_t dims1[2] = {M1->size1, M1->size2};
  herr_t status;
  
  fid    = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  status = H5LTmake_dataset(fid, "mat1", 2, dims1, H5T_NATIVE_DOUBLE, M1->data);
  assert(status != -1);
  
  hsize_t dims2[2] = {M2->size1, M2->size2};
  status = H5LTmake_dataset(fid, "mat2", 2, dims2, H5T_NATIVE_DOUBLE, M2->data);
  assert(status != -1);
  
  // Set some attributes:
  H5LTset_attribute_int(fid, "mat1", "attr1", &att1, (size_t) sizeof(int));
  
  status = H5Fclose(fid); assert(status != -1);
  assert(status != -1);
  
  return;
}

void hdf_read_two_matrix_main(gsl_matrix* M1, gsl_matrix* M2, int* attr, const char* fname) {
  hid_t  fid;
  herr_t status;
  
  outlog("Reading from %s...\n", fname);
  fid = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  
  H5LTread_dataset_double(fid, "mat1", M1->data);
  H5LTread_dataset_double(fid, "mat2", M2->data);
  H5LTget_attribute_int(fid, "mat1", "attr1", attr);
  
  outlog("done\n%s", "");
  status = H5Fclose(fid);
  assert(status != -1);

  return;
}

void hdf_read_matrix_main(gsl_matrix* M, const char* fname) {
  hid_t  fid;
  herr_t status;
  
  outlog("Reading from %s...\n",fname);
  fid    = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  H5LTread_dataset_double(fid, "term1", M->data);
  outlog("done\n%s","");
  status = H5Fclose(fid);
  assert(status != -1);

  return;
}

// Solves the linear system Ax = b for x.
// Assumes that x is already allocated.
void CholeskySolve_main(gsl_matrix* A, const gsl_vector* b, gsl_vector* x) {
  gsl_linalg_cholesky_decomp(A);
  gsl_linalg_cholesky_solve(A, b, x);
}

void LUSolve_main(gsl_matrix* A, const gsl_vector* b, gsl_vector* x) {
  int permutation_sign;
  gsl_permutation* permutation = gsl_permutation_alloc(b->size);

  gsl_linalg_LU_decomp(A, permutation, &permutation_sign);
  gsl_linalg_LU_solve(A, permutation, b, x);

  gsl_permutation_free(permutation);
}

int mkl_cblas_dgemm_main(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB,
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
#endif

int gsl2csr_main(gsl_matrix* D, mkl_spmat* &S) {
  const MKL_INT n = D->size1;
  const MKL_INT m = D->size2;
  S = (mkl_spmat*) malloc(sizeof(mkl_spmat));
  
  S->NZ = 0;
  S->n  = n;
  S->m  = m;
  S->Mi = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * (n+1), 64);
  MKL_INT    info; // Output status from MKL
  
  for (long i = 0; i < n*m; i++)
    if (D->data[i] != 0.0)
      S->NZ = S->NZ+1;
 
  MKL_INT  job[6] = {0, 0, 1, 2, S->NZ, 1};
  S->Mv  = (double*)  mklspace::mkl_malloc(sizeof(double)  * S->NZ, 64);
  S->Mj  = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * S->NZ, 64);

  mklspace::mkl_ddnscsr(job, &n, &m, D->data, &m, S->Mv, S->Mj, S->Mi, &info);

  return info;
}

int mkl_csradd_main(mkl_spmat* A, mkl_spmat* B, mkl_spmat* &C) {
  char   trans  = 'N';
  MKL_INT  req  = 1;
  double  beta  = 1.0;
  MKL_INT m     = A->m;
  MKL_INT n     = A->n;
  MKL_INT sort  = 0;
  MKL_INT nzmax = A->NZ + B->NZ;
  MKL_INT info;

  C = (mkl_spmat*) malloc(sizeof(mkl_spmat));
  
  C->n  = n;
  C->m  = m;
  C->Mi = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * (n+1), 64);
  
  mklspace::mkl_dcsradd(&trans, &req,  &sort, &m, &n, A->Mv, A->Mj, A->Mi,  &beta,
			B->Mv,  B->Mj, B->Mi, C->Mv,  C->Mj, C->Mi, &nzmax, &info);

  C->NZ = C->Mi[n]-1;
  C->Mv = (double*)  mklspace::mkl_malloc(sizeof(double)  * C->NZ, 64);
  C->Mj = (MKL_INT*) mklspace::mkl_malloc(sizeof(MKL_INT) * C->NZ, 64);

  req = 2;
  mklspace::mkl_dcsradd(&trans, &req,  &sort, &m, &n, A->Mv, A->Mj, A->Mi,  &beta,
			B->Mv,  B->Mj, B->Mi, C->Mv,  C->Mj, C->Mi, &nzmax, &info);

  return info;
}

int mkl_csrmul_main(mkl_spmat* A, mkl_spmat* B, mkl_spmat* &C) {
  char   trans  = 'N';
  MKL_INT  req  = 1;
  MKL_INT n     = A->n;
  MKL_INT m     = A->m;
  MKL_INT k     = B->n;
  MKL_INT sort  = 7;
  MKL_INT nzmax = A->NZ + B->NZ;
  MKL_INT info;
  
  C = (mkl_spmat*) malloc(sizeof(mkl_spmat));
  
  C->n  = k;
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


void compcol2csr31_main(mkl_spmat* S, gsl_spmatrix* C) {
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

void print_spmatrix(gsl_spmatrix* M) {
  for (int i = 0 ; i < M->size1 ; i++) {
    if (i==0)
      printf("0        1        2        3\n");
    for (int j = 0 ; j < M->size2 ; j++)
      printf("%f ",gsl_spmatrix_get(M, i, j));
    printf(" %i\n",i);
  }
      
  printf("\n(NZ=%i): \ndata     | i | p\n",M->nz);
  for (int i = 0 ; i < M->nz ; i++) {
    if (i <= M->size1) 
      printf("%f | %i | %i\n",M->data[i], M->i[i], M->p[i]);
    else
      printf("%f | %i | X\n",M->data[i], M->i[i]);
  }

  return;
}

/*
 * main function
 * supports fitting a dynamic topic model
 *
 */
int main(int argc, char* argv[]) {
    //   Various system-level options might need to be set...
    //   would be best put in a ./configure
    set_fpu(0x27F);
    //#pragma STDC FENV_ACCESS ON
    fesetround(FE_TONEAREST); // VERY IMPORTANT!!
    set_thread_affinity();
    gsl_ieee_env_setup();

	/* OLD WAY 
    //struct sigaction sa, osa;
    //sa.sa_flags = SA_ONSTACK | SA_RESTART | SA_SIGINFO;
    //sa.sa_sigaction = signal_handler;
    //sigaction(SIGFPE, &sa, &osa);
    omp_set_max_active_levels(2);
    omp_set_dynamic(1);
    mklspace::mkl_set_dynamic(0);
    mklspace::mkl_set_num_threads(FLAGS_threads);
	 END OLD WAY */

    omp_set_nested(true);
    omp_set_max_active_levels(2);
    omp_set_dynamic(0);
    omp_set_num_threads(FLAGS_threads);
    //mklspace::mkl_set_num_threads(FLAGS_threads);
    
    mklspace::mkl_set_interface_layer(MKL_INTERFACE_ILP64+MKL_INTERFACE_GNU);
    //mklspace::mkl_set_interface_layer(MKL_INTERFACE_ILP64);
#ifdef USETBB
    mklspace::mkl_set_threading_layer(MKL_THREADING_TBB);
#else
    mklspace::mkl_set_threading_layer(MKL_THREADING_INTEL);
#endif
    mklspace::mkl_set_dynamic(1);
    mklspace::mkl_set_num_threads(FLAGS_threads);

    google::ParseCommandLineFlags(&argc, &argv, 0);

#pragma omp parallel for num_threads(FLAGS_threads) schedule(static)
    for (int i = 0; i < FLAGS_threads; i++) {
	if (i ==0 ) {
    	   outlog("IOMP(num_procs)=%i\nIOMP(num_threads)=%i\nIOMP(max_threads)=%i\n",(int) omp_get_num_procs(), (int) omp_get_num_threads(), (int) omp_get_max_threads());
	}
	outlog("This is thread %i on core %i\n",(int) omp_get_thread_num(), sched_getcpu());
    }

    outlog("primitive int has size=%i chars of size=%i\n",sizeof(int),CHAR_BIT);

    // I know this seems outrageous, but one call to dgemm tends to work, maybe because of some 
    // kind of LD caching, but without the correctly linked interfaces (which will sometimes
    // not throw ld errors), the *second* call to blas will seg-fault. ie. Not the first.
    // Oh, and re-using the same variables caches the stack's function pointers, so it needs
    // new variables... FML
    outlog("Testing your BLAS library...%s","");
    gsl_matrix* D = gsl_matrix_alloc(2,2);
    gsl_matrix* E = gsl_matrix_alloc(2,2);
    gsl_matrix* F = gsl_matrix_alloc(2,2);

    for (int i = 0; i < 2; i++) 
      for (int j = 0; j < 2; j++) { 
	gsl_matrix_set(D, i, j, 0.1);
	gsl_matrix_set(E, i, j, 0.2);
      }
    
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, D, E, 0.0, F);

    gsl_matrix_free(D);
    gsl_matrix_free(E);
    gsl_matrix_free(F);

    gsl_matrix* M  = gsl_matrix_alloc(2,2);
    gsl_matrix* A  = gsl_matrix_alloc(2,2);
    gsl_matrix* B  = gsl_matrix_alloc(2,2);

    for (int i = 0; i < 2; i++) 
      for (int j = 0; j < 2; j++) {
	gsl_matrix_set(A, i, j, 0.1);
	gsl_matrix_set(B, i, j, 0.2);
      }
    
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, B, 0.0, M);

    gsl_matrix_free(M);
    gsl_matrix_free(B);
    gsl_matrix_free(A);

    outlog("success%s\n","");

#ifdef GTESTMODE
    if (FLAGS_model == "htest") {
      printf("In htest\n");
      gsl_matrix* S1 = gsl_matrix_alloc(5,5);
      char fname[BUFSIZ];
      sprintf(fname, "test.h5");

      gsl_matrix_set(S1,0,0,1.0);
      gsl_matrix_set(S1,1,1,2.0);
      gsl_matrix_set(S1,2,2,3.0);
      gsl_matrix_set(S1,3,3,4.0);
      gsl_matrix_set(S1,1,0,5.0);
      gsl_matrix_set(S1,0,1,6.0);
      gsl_matrix_set(S1,3,0,7.0);
      gsl_matrix_set(S1,0,3,8.0);

      hdf_write_matrix_main(S1, fname);

      gsl_matrix_free(S1);
      return 0;
    } 
    if (FLAGS_model == "stest") {
      int N=20;
      gsl_spmatrix* S1 = gsl_spmatrix_alloc(N,N);
      gsl_spmatrix* S2 = gsl_spmatrix_alloc(N,N);

      gsl_spmatrix_set(S1,0,0,1.0);
      gsl_spmatrix_set(S1,1,1,2.0);
      gsl_spmatrix_set(S1,2,2,3.0);
      gsl_spmatrix_set(S1,3,3,4.0);
      gsl_spmatrix_set(S1,1,0,5.0);
      gsl_spmatrix_set(S1,0,1,6.0);
      gsl_spmatrix_set(S1,3,0,7.0);
      gsl_spmatrix_set(S1,0,3,8.0);

      gsl_spmatrix_set(S2,0,0,1.0);
      gsl_spmatrix_set(S2,1,1,2.0);
      gsl_spmatrix_set(S2,2,2,3.0);
      gsl_spmatrix_set(S2,3,3,4.0);
      gsl_spmatrix_set(S2,1,0,5.0);
      gsl_spmatrix_set(S2,0,1,6.0);
      gsl_spmatrix_set(S2,3,0,7.0);
      gsl_spmatrix_set(S2,0,3,8.0);
      
      gsl_spmatrix* C1 = gsl_spmatrix_compcol(S1);
      gsl_spmatrix* C2 = gsl_spmatrix_compcol(S2);
      gsl_spmatrix* P  = gsl_spmatrix_alloc_nzmax(N,N,1,GSL_SPMATRIX_CCS);	
      
      printf("\nC1 = \n");
      print_spmatrix(C1);
      gsl_spblas_dgemm(1.0, C1, C2, P);
      printf("\nP = \n");
      print_spmatrix(P);
      
      gsl_spmatrix_free(S1);
      gsl_spmatrix_free(S2);
      gsl_spmatrix_free(P);
      gsl_spmatrix_free(C1);
      gsl_spmatrix_free(C2);
     
      return 1;
    }
    if (FLAGS_model == "ctest") {
      int N=4;
      gsl_spmatrix* S  = gsl_spmatrix_alloc(N,N);
      gsl_spmatrix* Sp = gsl_spmatrix_alloc(N,N);
      mkl_spmat*    Sm = (mkl_spmat*) malloc(sizeof(mkl_spmat));
      
      gsl_spmatrix_set(S,0,0,1.0);
      gsl_spmatrix_set(S,1,1,2.0);
      gsl_spmatrix_set(S,2,2,3.0);
      gsl_spmatrix_set(S,3,3,4.0);
      gsl_spmatrix_set(S,1,0,5.0);
      gsl_spmatrix_set(S,0,1,6.0);
      gsl_spmatrix_set(S,3,0,7.0);
      gsl_spmatrix_set(S,0,3,8.0);

      gsl_spmatrix* C = gsl_spmatrix_compcol(S);
      
      printf("S = \n");
      for (int i = 0 ; i < N ; i++) {
	for (int j = 0 ; j < N ; j++)
	  printf("%f ",gsl_spmatrix_get(S, i, j));
	printf(" %i\n",i);
      }

      printf("S data = \n data    | i | p\n");
      for (int i = 0 ; i < S->nz ; i++)
	printf("%f | %i | %i\n",S->data[i], S->i[i], S->p[i]);
      	
      printf("C data = \n data    | i | p\n");
      for (int i = 0 ; i < C->nz ; i++)
	printf("%f | %i | %i\n",C->data[i], C->i[i],C->p[i]);

      compcol2csr31_main(Sm,C);

      printf("Sm data = \n   Mv    | Mi | Mj\n");
      for (int i = 0 ; i < Sm->NZ ; i++) {
	if (i < (Sm->n)+1)
	  printf("%f | %i | %i\n",Sm->Mv[i], Sm->Mi[i],Sm->Mj[i]);
	else
	  printf("%f | X | %i\n",Sm->Mv[i],Sm->Mj[i]);
      }
      gsl_spmatrix_free(S);
      gsl_spmatrix_free(C);
      return 1;
    }
    else if (FLAGS_model == "ttest") {
      int n1 = 4;
      int n2 = 3;
      gsl_matrix* A  = gsl_matrix_calloc(n1,n2);
      gsl_matrix* B  = gsl_matrix_calloc(n2,n1);
      gsl_matrix* C  = gsl_matrix_calloc(n1,n1);
      mkl_spmat*  As = (mkl_spmat*) malloc(sizeof(mkl_spmat));
      mkl_spmat*  Bs = (mkl_spmat*) malloc(sizeof(mkl_spmat));
      mkl_spmat*  Cs;
      
      // A mostly sparse matrix:
      gsl_matrix_set(A, 0,         0,         1.0);
      gsl_matrix_set(A, 0,         1,         2.0);
      gsl_matrix_set(A, n1,        n2,         3.0);
      gsl_matrix_set(B, 0,         0,         1.0);
      gsl_matrix_set(B, 0,         1,         2.0);
      gsl_matrix_set(B, 1,         1,         2.0);
      gsl_matrix_set(B, n2,         n1,         3.0);

      gsl2csr_main(A,As);
      gsl2csr_main(B,Bs);
      
      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,B,0.0,C);
      printf("GSL Result=\n");
      for (int i = 0 ; i < n1; i++) {
	for (int j = 0 ; j < n1; j++)
	  printf("%f ",gsl_matrix_get(C, i, j));
	printf(" %i\n",i);
      }

      //      mkl_csradd_main(As,Bs,Cs);
      mkl_csrmul_main(As,Bs,Cs);
      printf("NZ C=%i\n",Cs->NZ);

      for (int i = 0; i<Cs->NZ; i++)
	printf("Cv[%i]=%f\n",i,Cs->Mv[i]);
      for (int i = 0; i<Cs->NZ; i++)
	printf("Ci[%i]=%i\n",i,Cs->Mi[i]);
      for (int i = 0; i<Cs->NZ; i++)
	printf("Cj[%i]=%i\n",i,Cs->Mj[i]);

      gsl_matrix_free(A);
      gsl_matrix_free(B);
      
      return 0;
    }
    else if (FLAGS_model == "sptest") {
      int N  = 10000;
      int n1 = N;
      int n2 = N;
      gsl_matrix* A = gsl_matrix_calloc(n1, n2);
      gsl_matrix* B = gsl_matrix_calloc(n1, n2);
      gsl_matrix* M = gsl_matrix_calloc(n2, n2);
      
      // A mostly sparse matrix:
      outlog("SP 1%s\n","");
      gsl_matrix_set(A, 0,         (int) N/2, 1.0);
      gsl_matrix_set(A, (int) N/3, 1,         2.0);
      gsl_matrix_set(A, 1,         (int) N/4, 3.0);
      gsl_matrix_set(B, 0,         (int) N/5, 4.0);
      gsl_matrix_set(B, (int) N/2, 0,         5.0);
      gsl_matrix_set(B, (int) N/3, 1,         6.0);
      gsl_matrix_set(B, 0,         (int) N/4, 7.0);
      outlog("SP 2%s\n","");

      gsl_matrix* Ac = gsl_matrix_calloc(n1,n2);
      gsl_matrix* Bc = gsl_matrix_calloc(n1,n2);
      gsl_matrix* Mc = gsl_matrix_calloc(n2,n2);
      gsl_matrix_memcpy(Ac,A);
      gsl_matrix_memcpy(Bc,B);

      outlog("SP 3%s\n","");
      clock_t gslt,mklt,spmklt,startt = clock();
      //      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, B, 0.0, M);
      gslt = clock()-startt;
      outlog("SP 4%s\n","");

      gsl_matrix_memcpy(Mc,M);
      assert(gsl_matrix_equal(Ac,A));
      assert(gsl_matrix_equal(Bc,B));
      outlog("SP 5%s\n","");

      M  = gsl_matrix_calloc(n2,n2);
      startt = clock();
      //      mkl_cblas_dgemm_main(CblasTrans, CblasNoTrans, 1.0, A, B, 0.0, M);
      mklt = clock()-startt;
      assert(gsl_matrix_equal(Ac,A));
      assert(gsl_matrix_equal(Bc,B));
      assert(gsl_matrix_equal(Mc,M));
      outlog("SP 6 *********%s\n","");

      M  = gsl_matrix_calloc(n2,n2);
      startt = clock();
      mkl_spblas_dgemm_main(A,B,M);
      outlog("SP 7%s\n","");

      spmklt = clock()-startt;
      assert(gsl_matrix_equal(Ac,A));
      assert(gsl_matrix_equal(Bc,B));
      assert(gsl_matrix_equal(Mc,M));

      outlog("**   GSL: %d\n**   MKL: %d\n** SPMKL: %d\n",gslt,mklt,spmklt);
      //      outlog("Free 1%s\n","");
      gsl_matrix_free(A);
      /* outlog("Free 2%s\n",""); */
      gsl_matrix_free(B);
      /* outlog("Free 3%s\n",""); */
      gsl_matrix_free(M);
      
      return 0;
    }
    else if (FLAGS_model == "choltest") {
      int         N = 9500;
      gsl_matrix *B = gsl_matrix_alloc(N,N);
      gsl_matrix *C = gsl_matrix_alloc(N,N);
      gsl_matrix *A = gsl_matrix_alloc(N,N);
      gsl_vector *b = gsl_vector_alloc(N);
      gsl_vector *x = gsl_vector_alloc(N);

      outlog("Setting%s\n","");
      for (int i = 0 ; i < N; i++)  {
	gsl_vector_set(b, i , (double) (i / 2));

	for (int j = 0 ; j < N; j++) {
	  gsl_matrix_set(B, i, j, i+j / 4);
	  gsl_matrix_set(C, i, j, i+j / 4);
	}
      }

      outlog("Mult%s\n","");
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, C, B, 0.0, A);
      outlog("Solve%s\n","");
      MKLSolve_lapacke_single_main(A,b,x);

      gsl_vector_free(x);
      gsl_vector_free(b);
      gsl_matrix_free(A);
      gsl_matrix_free(C);
      gsl_matrix_free(B);
      
      return 0;
    }
    else if (FLAGS_model == "dgemmtest") {
      int N = 5500;
      gsl_matrix *A = gsl_matrix_alloc(N,N);
      gsl_matrix *B = gsl_matrix_alloc(N,N);
      gsl_matrix *C = gsl_matrix_alloc(N,N);

      for (int i = 0 ; i < N; i++) 
	for (int j = 0 ; j < N; j++) 
	  gsl_matrix_set(A, i, j, i+j / 3);

      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, B, 0.0, C);
            
      gsl_matrix_free(A);
      gsl_matrix_free(B);
      gsl_matrix_free(C);
      return 0;
    }
    else if (FLAGS_model == "hdftest2") {
      gsl_matrix *A = gsl_matrix_alloc(4000,4000);
      
      char name[400];
      sprintf(name, "./dset.h5");
      hdf_read_matrix_main(A, name);
      
      gsl_matrix_free(A);
      return 0;
    }
    else if (FLAGS_model == "hdftest") {
      gsl_matrix *A = gsl_matrix_alloc(4000,4000);
      gsl_matrix *B = gsl_matrix_alloc(2000,2000);
      
      for (int i = 0 ; i < 1000; i++) 
	for (int j = 0 ; j < 1000; j++) 
	  gsl_matrix_set(A, i, j, i+j / 3);

      for (int i = 0 ; i < 1000; i++) 
	for (int j = 0 ; j < 1000; j++) 
	  gsl_matrix_set(A, i, j, i+j / 2);

      outlog("writing%s...","");
      char name[400];
      sprintf(name, "./double_dset.h5");  
      int attr1 = 7;
    
      hdf_write_two_matrix_main(A,B,attr1,name);
      outlog("done\n%s","");

      gsl_matrix_free(A);
      gsl_matrix_free(B);
      
      gsl_matrix *A2 = gsl_matrix_alloc(4000,4000);
      gsl_matrix *B2 = gsl_matrix_alloc(2000,2000);
      int attr2;

      hdf_read_two_matrix_main(A2,B2,&attr2,name);
      outlog("Read A1 dims = %i x %i\nB1 dims = %i x %i\nattr = %i\n",
	     A2->size1, A2->size2, B2->size1, B2->size2, attr2);

      return 0;
    }
    else if (FLAGS_model == "mkltest") {
      int NRHS = 1;
      gsl_matrix *A = gsl_matrix_alloc(4,4);
      gsl_vector *b = gsl_vector_alloc(4);
      gsl_vector *x = gsl_vector_alloc(4);
      
      gsl_matrix_set(A,0,0, 1.0);
      gsl_matrix_set(A,0,1, 1.0);
      gsl_matrix_set(A,0,2, 2.0);
      gsl_matrix_set(A,0,3,-5.0);
      gsl_matrix_set(A,1,0, 2.0);
      gsl_matrix_set(A,1,1, 5.0);
      gsl_matrix_set(A,1,2,-1.0);
      gsl_matrix_set(A,1,3,-9.0);
      gsl_matrix_set(A,2,0, 2.0);
      gsl_matrix_set(A,2,1, 1.0);
      gsl_matrix_set(A,2,2,-1.0);
      gsl_matrix_set(A,2,3, 3.0);
      gsl_matrix_set(A,3,0, 1.0);
      gsl_matrix_set(A,3,1, 3.0);
      gsl_matrix_set(A,3,2, 2.0);
      gsl_matrix_set(A,3,3, 7.0); 
      
      gsl_vector_set(b, 0, 3.0);
      gsl_vector_set(b, 1, -3.0);
      gsl_vector_set(b, 2, -11.0);
      gsl_vector_set(b, 3, -5.0);
      
      int ret = MKLSolve_lapacke_single_main(A,b,x);
      
      printf("MKL Resulting x=\n");
      for (int i = 0 ; i < 4; i++)
	printf("%f\n",gsl_vector_get(x, i));

      printf("Resulting b=\n");
      for (int i = 0 ; i < 4; i++)
	printf("%f\n",gsl_vector_get(b, i));

      printf("Resulting A=\n");
      for (int i = 0 ; i < 4; i++) {
	for (int j = 0 ; j < 4; j++)
	  printf("%f ",gsl_matrix_get(A, i, j));
	printf("\n");
      }

      LUSolve_main(A, b, x);
      printf("LU\nResulting x=\n");
      for (int i = 0 ; i < 4; i++)
	printf("%f\n",gsl_vector_get(x, i));

      printf("Resulting b=\n");
      for (int i = 0 ; i < 4; i++)
	printf("%f\n",gsl_vector_get(b, i));

      printf("Resulting A=\n");
      for (int i = 0 ; i < 4; i++) {
	for (int j = 0 ; j < 4; j++)
	  printf("%f ",gsl_matrix_get(A, i, j));
	printf("\n");
      }

      CholeskySolve_main(A, b, x);
      printf("Cholesky\nResulting x=\n");
      for (int i = 0 ; i < 4; i++)
	printf("%f\n",gsl_vector_get(x, i));

      printf("Resulting b=\n");
      for (int i = 0 ; i < 4; i++)
	printf("%f\n",gsl_vector_get(b, i));

      printf("Resulting A=\n");
      for (int i = 0 ; i < 4; i++) {
	for (int j = 0 ; j < 4; j++)
	  printf("%f ",gsl_matrix_get(A, i, j));
	printf("\n");
      }


      gsl_vector_free(x); 
      gsl_vector_free(b);
      gsl_matrix_free(A);
      
      return 0;
    }
#endif

    if (FLAGS_model == "rdim") {
      outlog("Checking if rDIM meta-data is well-formed...%s\n","");
      corpus_seq_t* c = read_corpus_seq(FLAGS_corpus_prefix.c_str(), FLAGS_ntopics, FLAGS_metafields);
      free_corpus_seq(c);

      // If we got a held-out corpus, check it too...
      if (strcmp("",FLAGS_heldout_corpus_prefix.c_str()) != 0) {
	corpus_seq_t* h = read_corpus_seq(FLAGS_heldout_corpus_prefix.c_str(), FLAGS_ntopics, FLAGS_metafields);
	free_corpus_seq(h);
      }

      outlog("Sequencing, meta-data and LDAC data all look good!%s\n", "");
    }

    // mode for spitting out document sums
    if (FLAGS_mode == "sums") {
      corpus_seq_t* c = read_corpus_seq(FLAGS_corpus_prefix.c_str(), FLAGS_ntopics, FLAGS_metafields);
      outlog("Read corpus %s\n", FLAGS_corpus_prefix.c_str());

      for (int t = 0; t < c->len; t++) {
	int sum = 0;
	for (int d = 0; d < c->corpus[t]->ndocs; d++)
	  sum += c->corpus[t]->doc[d]->total;

	printf("%d\n\n", sum);
      }
    }

    // mode for fitting a dynamic topic model
    else if (FLAGS_mode == "fit")
      fit_dtm(0, FLAGS_heldout_time - 1);

    // mode for analyzing documents through time according to a DTM
    // AG: This mode is undocumented, not sure what it does...
    else if (FLAGS_mode == "time") {
      // load corpus and model based on information from params
      corpus_seq_t* data = read_corpus_seq(FLAGS_heldout_corpus_prefix.c_str(), 
					   FLAGS_ntopics, FLAGS_metafields);
      lda_seq*     model = read_lda_seq(FLAGS_lda_model_prefix.c_str(), data);

      // initialize the table (D X OFFSETS)
      double**     table = (double**) malloc(sizeof(double*) * data->len);
      
      for (int t = 0; t < data->len; t++) {
	table[t] = (double*) malloc(sizeof(double) * data->corpus[t]->ndocs);

	for (int d = 0; d < data->corpus[t]->ndocs; d++)
	  table[t][d] = -1;  // this should be NAN
      }

      // set the LDA model to be populated
      lda_post post;
      lda* lda_model = new_lda_model(model->ntopics, model->nterms);
      int max_nterms = compute_max_nterms(data);
      post.phi       = gsl_matrix_calloc(max_nterms, model->ntopics);
      post.log_phi   = gsl_matrix_calloc(max_nterms, model->ntopics);
      post.gamma     = gsl_vector_calloc(model->ntopics);
      post.lhood     = gsl_vector_calloc(model->ntopics);
      post.model     = lda_model;

      // compute likelihoods for each model
      for (int t = 0; t < data->len; t++) {
	make_lda_from_seq_slice(lda_model, model, t);

	for (int d = 0; d < data->corpus[t]->ndocs; d++) {
	  post.doc = data->corpus[t]->doc[d];
	  double likelihood = fit_lda_post(d, t, &post, model, NULL, NULL, NULL, NULL);
	  table[t][d] = post.doc->log_likelihood;
	}
      }

      char tmp_string[400];
      sprintf(tmp_string, "%s-heldout_post.dat", FLAGS_outname.c_str());
      FILE* post_file = fopen(tmp_string, "w");
      
      for (int t=0; t < data->len; ++t) {
	if (data->corpus[t]->ndocs >= 0)
	  fprintf(post_file, "%f", table[t][0]);
	
	for (int d = 1; d < data->corpus[t]->ndocs; ++d)
	  fprintf(post_file, ",%f", table[t][d]);
	
	fprintf(post_file, "\n");
      }
    }
    else
      outlog("Invalid mode '%s'\n", FLAGS_mode.c_str());
    
#ifdef CUDA
    cusparseDestroy(hndl);
#endif
    
    return(0);
}
