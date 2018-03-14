//           Aaron Gerow (gerow@uchicago.edu)
//
#include <error.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <assert.h>
#include <omp.h>
#define GOMP
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <fenv.h>
#include "gsl/gsl_ieee_utils.h"
#include <gsl/gsl_randist.h>

#include "hdf5.h"
#include "hdf5_hl.h"
#include <math.h>
#include "gflags.h"

namespace mklspace
{
#include "mkl_lapacke.h"
#include "mkl.h"
#include "mkl_types.h"
#include "mkl_spblas.h"
#include "mkl_service.h"
#include "mkl_pardiso.h"
}

DEFINE_string(wd,"","Working directory.");
DEFINE_string(Aloc,"","Location of A in CSR3-HDF5 format.");
#define LARGE_TEST_SIZE 10000
#define NRHS 32

// SSE 4.5 enabled CPUs:
void set_fpu (unsigned int mode) {
  asm("fldcw %0" : : "m" (*&mode));
}

int hdf_read_csr3(double* &a, MKL_INT* &ia, MKL_INT* &ja, MKL_INT &NZ, MKL_INT &n1, const char* fname) {
  hid_t  fid;
  herr_t status;
  
  fprintf(stderr, "INFO: (in pardiso_solver): Reading CSR3 matrix...\n");
  fid    = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
  H5LTget_attribute_long_long(fid, "CSR3_a", "n1", &n1);
  H5LTget_attribute_long_long(fid, "CSR3_a", "NZ", &NZ);

  if (n1 < 1 || NZ < 1) {
    fprintf(stderr,"ERROR (in pardiso_solve): Bad CSR3 form (n+1=%i, NZ=%i): exiting...\n", fname, n1,NZ);
    _exit(24);
  }
  
  ia = (MKL_INT*) malloc(sizeof(MKL_INT) * n1);
  ja = (MKL_INT*) malloc(sizeof(MKL_INT) * NZ);
  a  = (double*)  malloc(sizeof(double)  * NZ); 

  status = H5LTread_dataset_double(   fid, "CSR3_a",   a);
  status = H5LTread_dataset(fid, "CSR3_ia", H5T_NATIVE_LLONG, ia);
  status = H5LTread_dataset(fid, "CSR3_ja", H5T_NATIVE_LLONG, ja);
  status = H5Fclose(fid);

  if (status == -1) {
    fprintf(stderr,"(ERROR (in pardiso_solver): HDF5 returned -1 on reading %s: %i\n", fname, (int) status);
    _exit(22);
  }
  else {
    fprintf(stderr,"done reading CSR3 (status=%i)\n",(int) status);
  }

  for (long i = 0 ; i < NZ ; i++) {
    if (ja[i] < 0)
      printf("PARDISO ERROR S1: ja[%i]=%i (v=%f)\n",i,ja[i],a[i]);
  }


  return (int) status;
}

int pardiso_solve(gsl_matrix* gsl_A, gsl_matrix* gsl_b, gsl_matrix* gsl_x) {
  MKL_INT 
    n  = (MKL_INT) gsl_A->size1,
    m  = (MKL_INT) gsl_A->size2,
    NZ = 0,     
    mtype = -2,          /* Real symmetric matrix */ // -2 is symmetric indeffinite
    nrhs = gsl_b->size2, /* Number of right hand sides. */
    iparm[64],
    maxfct, mnum, phase, error, msglvl, idum,
    c = 0, jc=0, col=0;
  double v = 0, ddum;
  void *pt[64];

  for (MKL_INT i=0; i<gsl_A->size1; i++) 
    for (MKL_INT j=0; j<i; j++) 
      gsl_matrix_set(gsl_A,i,j,0.0); 

  for (MKL_INT i=0;i<n*n;i++)
    if (gsl_A->data[i] != 0.0) 
      NZ++;

  MKL_INT
    *ia = (MKL_INT*) malloc(sizeof(MKL_INT) * (n+1)),
    *ja = (MKL_INT*) malloc(sizeof(MKL_INT) * NZ);  
  double 
    *b  = (double*) malloc(sizeof(double) * n * nrhs), 
    *x  = (double*) malloc(sizeof(double) * n * nrhs),
    *a  = (double*) malloc(sizeof(double) * NZ); 

  c=0;
  for (MKL_INT j=0; j<nrhs; j++)
    for (MKL_INT i = 0; i < n; i++ )
      b[c++] = gsl_matrix_get(gsl_b,i,j);
  c=0;

    bool rowset;
  for (MKL_INT i=0; i<n; i++) {
    rowset = false;
    for (MKL_INT j=i; j<m; j++) { // Note this is only the upper triangle
      v = gsl_matrix_get(gsl_A,i,j);

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

  for (int i = 0; i < 64; i++) {
    pt[i] = 0;
    iparm[i] = 0;
  }

  iparm[0]  = 1;  /* No solver default */
  iparm[1]  = 2;  /* Fill-in reordering from METIS */ // 3 would be an MP version...
  iparm[3]  = 0;  /* No iterative-direct algorithm */
  iparm[4]  = 0;  /* No user fill-in reducing permutation */
  iparm[5]  = 0;  /* 0 Write solution into x */
  iparm[7]  = 2;  /* Max numbers of iterative refinement steps */
  iparm[9]  = 13; /* Perturb the pivot elements with 1E-13 */
  iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
  iparm[11] = 0;  /* Solve system A X = B (as oppose to A^T X = B*/
  iparm[12] = 0;  /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
  iparm[13] = 0;  /* Output: Number of perturbed pivots */
  iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */ // Writes to screen
  iparm[18] = -1; /* Output: Mflops for LU factorization */  // writes to screen
  iparm[19] = 0;  /* Output: Numbers of CG Iterations */
  iparm[20] = 0;  /* Try one for very indefinite matrices */
  iparm[23] = 0;  /* 1 can improve factoring time when Threads >=~ 8 */ // Not sure on AMD or libgomp systems // AG
  iparm[26] = 0;  /* Check the CSR format on input */ // 0 Does not check
  iparm[27] = 0;  /* 0 is double precision, 1 is single precision */
  iparm[34] = 1;  /* 0 is 1-based indexing, 1 is 0-based */
  iparm[59] = 0;  /* Out-of-core processing disabled */ // 2 will use disk for permutations and pivot tabes.
  maxfct    = 1;  /* Maximum number of numerical factorizations. */
  mnum      = 1;  /* Which factorization to use. */
  msglvl    = 1;  /* Print statistical information in file */
  error     = 0;  /* Initialize error flag */

  phase = 11;
  mklspace::PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		     &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
  if ( error != 0 )    {
    printf ("\nERROR during symbolic factorization: %d", error);
    exit (1);
  }
  printf ("\nReordering completed ... "); 
  printf ("\nNumber of nonzeros in factors = %d", iparm[17]); 
  printf ("\nNumber of factorization MFLOPS = %d", iparm[18]); 

  phase = 22;
  mklspace::PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		     &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
  if ( error != 0 )    {
    printf ("\nERROR during numerical factorization: %d", error);
    exit (2);
  }
  iparm[7] = 2;         // Max numbers of iterative refinement steps. 

  phase = 23;
  mklspace::PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		     &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
  if ( error != 0 )    {
    printf ("\nERROR during solution: %d", error);
    exit (3);
  }

  c=0;
  for (MKL_INT j=0; j<nrhs; j++)
    for (MKL_INT i=0; i<n;i++)
      gsl_matrix_set(gsl_x,i,j,x[c++]);

  phase = -1;           /* Release internal memory. */
  mklspace::PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		     &n, &ddum, ia, ja, &idum, &nrhs,
		     iparm, &msglvl, &ddum, &ddum, &error);

  printf("Done with PARDISO objects, freeing...\n");

  free(ia);
  free(ja);
  free(a);
  free(b);
  free(x);
}

int pardiso_solve_csr3(gsl_matrix* gsl_b, gsl_matrix* gsl_x, const char* fname) {
  double *a = NULL;
  MKL_INT *ja = NULL, *ia = NULL;
  MKL_INT NZ,n1;
  hdf_read_csr3(a,ia,ja,NZ,n1,fname);
  MKL_INT n = n1-1;

  MKL_INT 
    mtype = -2,          /* Real symmetric matrix */ // -2 is symmetric indeffinite
    nrhs = gsl_b->size2, /* Number of right hand sides. */
    iparm[64],
    maxfct, mnum, phase, error, msglvl, idum,
    c=0;
  double ddum;
  void *pt[64];

  double 
    *b  = (double*) malloc(sizeof(double) * n * nrhs), 
    *x  = (double*) malloc(sizeof(double) * n * nrhs);

  c=0;
  for (MKL_INT j=0; j<nrhs; j++)
    for (MKL_INT i = 0; i < n; i++ )
      b[c++] = gsl_matrix_get(gsl_b,i,j);
  c=0;
 
  for (int i = 0; i < 64; i++ ){
    pt[i] = 0;
    iparm[i] = 0;
  }

  iparm[0]  = 1;  /* No solver default */
  iparm[1]  = 2;  /* Fill-in reordering from METIS */ // 3 would be an MP version...
  iparm[3]  = 0;  /* No iterative-direct algorithm */
  iparm[4]  = 0;  /* No user fill-in reducing permutation */
  iparm[5]  = 0;  /* 0 Write solution into x */
  iparm[7]  = 2;  /* Max numbers of iterative refinement steps */
  iparm[9]  = 13; /* Perturb the pivot elements with 1E-13 */
  iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
  iparm[11] = 0;  /* Solve system A X = B (as oppose to A^T X = B*/
  iparm[12] = 0;  /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
  iparm[13] = 0;  /* Output: Number of perturbed pivots */
  iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */ // Writes to screen
  iparm[18] = -1; /* Output: Mflops for LU factorization */  // writes to screen
  iparm[19] = 0;  /* Output: Numbers of CG Iterations */
  iparm[20] = 0;  /* Try one for very indefinite matrices */
  iparm[23] = 0;  /* 1 can improve factoring time when Threads >=~ 8 */ 
  iparm[26] = 1;  /* Check the CSR format on input */ // 0 Does not check
  iparm[27] = 0;  /* 0 is double precision, 1 is single precision */
  iparm[34] = 1;  /* PARDISO use C-style indexing for ia and ja arrays */
  iparm[59] = 0;  /* Out-of-core processing disabled */ // 2 will use disk for permutations and pivot tabes.
  maxfct    = 1;  /* Maximum number of numerical factorizations. */
  mnum      = 1;  /* Which factorization to use. */
  msglvl    = 1;  /* Print statistical information in file */
  error     = 0;  /* Initialize error flag */

  phase = 11;
  mklspace::PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		     &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
  if ( error != 0 )    {
    printf ("\nERROR during symbolic factorization: %d", error);
    exit (1);
  }
  printf ("\nReordering completed ... "); 
  printf ("\nNumber of nonzeros in factors = %d", iparm[17]); 
  printf ("\nNumber of factorization MFLOPS = %d", iparm[18]); 

  phase = 22;
  mklspace::PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		     &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
  if ( error != 0 )    {
    printf ("\nERROR during numerical factorization: %d", error);
    exit (2);
  }
  iparm[7] = 2;         // Max numbers of iterative refinement steps. 

  phase = 23;
  mklspace::PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		     &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
  if ( error != 0 )    {
    printf ("\nERROR during solution: %d", error);
    exit (3);
  }

  c=0;
  for (MKL_INT j=0; j<nrhs; j++)
    for (MKL_INT i=0; i<n;i++)
      gsl_matrix_set(gsl_x,i,j,x[c++]);

  phase = -1;           /* Release internal memory. */
  mklspace::PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		     &n, &ddum, ia, ja, &idum, &nrhs,
		     iparm, &msglvl, &ddum, &ddum, &error);

  printf("Done with PARDISO objects, freeing...\n");

  free(ia);
  free(ja);
  free(a);
  free(b);
  free(x);
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


gsl_matrix* hdf_read_matrix(const char* fname) {
  gsl_matrix*   ret;
  hid_t         fid;
  herr_t        status;
  unsigned long s1, s2;
  
  fid    = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
  H5LTget_attribute_ulong(fid, "term1", "size1", &s1);
  H5LTget_attribute_ulong(fid, "term1", "size2", &s2);

  if (s1 < 1 || s2 < 1) {
    fprintf(stderr,"ERROR (in pardiso_solve): Got matrix->size1 < 1 for %s (%i x %i), exiting...\n", fname, s1,s2);
    _exit(24);
  }

  ret = gsl_matrix_alloc(s1, s2);

  status = H5LTread_dataset_double(fid, "term1", ret->data);
  status = H5Fclose(fid);

  if (status == -1) {
    fprintf(stderr,"(ERROR (in pardiso_solve): HDF5 returned -1 on reading %s: %i\n", fname, (int) status);
    _exit(22);
  }

  return ret;
}

void run_mrhs_test() {
  int         N = LARGE_TEST_SIZE;
  int      nrhs = NRHS;
  gsl_matrix *B = gsl_matrix_calloc(N,N); 
  gsl_matrix *A = gsl_matrix_calloc(N,N);
  gsl_matrix *b = gsl_matrix_calloc(N,nrhs);
  gsl_matrix *x = gsl_matrix_calloc(N,nrhs);
  gsl_rng    *r = gsl_rng_alloc(gsl_rng_taus);
  time_t startt;
  double val;

  gsl_rng_set(r, 1); // Seed
  int j, N_n=N/4;

  printf("Running larger test %i x %i with %i rhs.\n",N,N,nrhs);

  for (int rhs=0; rhs<nrhs; rhs++) {
    for (int i = 0 ; i < N; i++)  {
      val = gsl_ran_gaussian(r, 0.1); // Second arg is \sigma in N(0,\sigma)
      gsl_matrix_set(b, i, rhs,(val*i));
      gsl_matrix_set(B, i, i, (double) (1+i+rhs / (double) N)); 
    }
  }

  // Add noise to N/4 values in each row
  for (int i=0; i<N; i++) {
    do {
      val = gsl_ran_gaussian(r, 0.1); // Second arg is \sigma in N(0,\sigma)
      j   = gsl_ran_flat(r,0,N_n);
      gsl_matrix_set(B,i,j,val);
    } while (N_n-- > 0);
  }

  gsl_matrix_scale(B,1.0/N);
  gsl_matrix_scale(b,1.0/N);

  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, B, B, 0.0, A); 

  printf("\nInitial b[0:4::0:4] = \n");
  for (int i = 0 ; i < 4; i++) { 
  for (int j = 0 ; j < 4; j++)
    printf("%f ",gsl_matrix_get(b, i, j));
  printf("\n");
  }

  printf("\nInitial A[0:4::0:4]=\n"); 
  for (int i = 0 ; i < 4; i++) { 
  for (int j = 0 ; j < 4; j++) 
    printf("%f ",gsl_matrix_get(A, i, j)); 
  printf("\n"); 
  } 

  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, B, B, 0.0, A); 
  startt=clock();
  pardiso_solve(A,b,x); // Clobbers A ?and maybe b?

  printf("\nMKL Result x[0:4::0:4] (%i seconds) = \n",clock()-startt);
  for (int i = 0 ; i < 4; i++) { 
  for (int j = 0 ; j < 4; j++)
    printf("%f ",gsl_matrix_get(x, i, j));
  printf("\n");
  }

  gsl_matrix_free(x); 
  gsl_matrix_free(b);
  gsl_matrix_free(A);
  gsl_matrix_free(B);
}

void run_aps_test() {
  char fname[BUFSIZ];
  sprintf(fname, "/scratch/midway/gerow/outputs-project/aps-from_iter0/A.h5");
  gsl_matrix *A = hdf_read_matrix(fname);
  sprintf(fname, "/scratch/midway/gerow/outputs-project/aps-from_iter0/b.h5");
  gsl_matrix *b = hdf_read_matrix(fname);

  int         N = b->size1;
  int      nrhs = b->size2;
  gsl_matrix *x = gsl_matrix_calloc(N,nrhs);

  printf("\nInitial b[0:4::0:4] = \n");
  for (int i = 0 ; i < 4; i++) { 
    for (int j = 0 ; j < 4; j++)
      printf("%f ",gsl_matrix_get(b, i, j));
    printf("\n");
  }

  printf("\nInitial A[0:4::0:4]=\n"); 
  for (int i = 0 ; i < 4; i++) { 
    for (int j = 0 ; j < 4; j++) 
      printf("%f ",gsl_matrix_get(A, i, j)); 
    printf("\n"); 
  } 

  pardiso_solve(A,b,x); // Clobbers A ?and maybe b?

  printf("\nMKL Result x[0:4::0:4] = \n");
  for (int i = 0 ; i < 4; i++) { 
    for (int j = 0 ; j < 4; j++)
      printf("%f ",gsl_matrix_get(x, i, j));
    printf("\n");
  }

  gsl_matrix_free(x); 
  gsl_matrix_free(b);
  gsl_matrix_free(A);
}

int main(int argc, char* argv[]) {
  set_fpu(0x27F);
  //#pragma STDC FENV_ACCESS ON
  fesetround(FE_TONEAREST);

  google::ParseCommandLineFlags(&argc, &argv, 0);
  fprintf(stderr, "*** This is pardiso_solve ***\n");
  fprintf(stderr, "GOMP(num_procs)=%i\nGOMP(max_threads)=%i\n",(int) omp_get_num_procs(), (int) omp_get_max_threads());
  fprintf(stderr, "wd = %s\n",FLAGS_wd.c_str());
  fprintf(stderr, "G_PARDISO_SOLVE: Reading HDF5 b...\n");

  if (false) {
    run_aps_test();
    fprintf(stderr, "DONE WITH TEST!!\n");
    sleep(30);
    _exit(-1);
  }
  
  if (false) {
    run_mrhs_test();
    fprintf(stderr, "DONE WITH TEST!!\n");
    sleep(40);
    _exit(-1);
  }

  char fname[BUFSIZ];

  sprintf(fname, "%s/b.h5",FLAGS_wd.c_str());
  gsl_matrix* b = hdf_read_matrix(fname); // b doesn't get clobbered by Solves
  fprintf(stderr, "done\nAllocating x...");
  gsl_matrix* x = gsl_matrix_alloc(b->size1,b->size2);

  fprintf(stderr, "done\nG_PARDISO_SOLVE: Expecting HDF5 version of A in CSR3 representation (%s)...\n", FLAGS_Aloc.c_str());
  fprintf(stderr, "Solving with MKL::PARDISO x...\n");
  
  sprintf(fname, "%s",FLAGS_Aloc.c_str());
  pardiso_solve_csr3(b,x,fname);

  fprintf(stderr, "\nHDF MKL::PARDISO Result x[0:4::0:4] = \n");
  for (int i = 0 ; i < 4; i++) { 
    for (int j = 0 ; j < 4; j++)
      fprintf(stderr, "%f ",gsl_matrix_get(x, i, j));
    fprintf(stderr, "\n");
  }

  fprintf(stderr, "\nWriting x...\n");
  sprintf(fname, "%s/x.h5",FLAGS_wd.c_str());
  hdf_write_matrix(x,fname);
  
  printf("Freeing A, b and x\n");
  gsl_matrix_free(x); 
  gsl_matrix_free(b);

  return 0;
}
