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


#ifndef LDASEQ_H
#define LDASEQ_H
#include <sys/stat.h>
#include <sys/types.h>
#include <float.h>
#include "gsl-wrappers.h"
#include "lda.h"

#define LDA_SEQ_EM_THRESH 1e-4
#define SAVE_LAG 10

/*
 * an lda sequence is a collection of simplex sequences for K topics
 * and an alpha vector
 *
 */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

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

#ifdef SPARSE
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>
#endif

#include <stdlib.h>
#include <assert.h>
#include "param.h"
#include "ss-lm.h"
#include "data.h"
#include "lda.h"
#define LDA_SEQ_EM_THRESHOLD 1e-5;

// lda sequence variational posterior distribution
typedef struct mkl_spmat {
  double*  Mv;
  MKL_INT* Mi;
  MKL_INT* Mj;
  MKL_INT  NZ;
  MKL_INT  n;
  MKL_INT  m;
} mkl_spmat;


// === allocation and initialization ===
inf_var* inf_var_alloc(int number_topics, corpus_seq_t* corpus_seq);
void inf_var_free(inf_var* ptr);
void mkl_spmat_free(mkl_spmat* ptr);

// initialize lda sequence from lda model topics
void init_lda_seq_from_ss(lda_seq* model, double      topic_chain_variance, double topic_obs_variance, 
			  double   alpha, gsl_matrix* init_suffstats);

// === fitting ===

// infer a corpus with an lda-seq
double update_inf_var(lda_seq*     seq, const corpus_seq_t* data, 
		      gsl_matrix** phi,              size_t t,   const char* root);
double update_inf_var_multiple(lda_seq*     seq, const corpus_seq_t* data, 
			       gsl_matrix** phi,              size_t t,   const char* root);
void   update_inf_reg(lda_seq*     seq, const corpus_seq_t* data, 
		      gsl_matrix** phi,              size_t t, const char* root);
double lda_seq_infer(lda_seq*     model,     const corpus_seq_t* data, 
		     gsl_matrix** suffstats,         gsl_matrix* gammas,
                     gsl_matrix*  lhoods,                    int iter, 
		     const char* file_root,           const bool heldout);

// fit lda sequence from sufficient statistics

double fit_lda_seq(lda_seq* model, const corpus_seq_t* data, const corpus_seq_t* heldout, const char* file_root);
void   update_lda_seq_ss(int time, const doc_t* doc, const lda_post* post, gsl_matrix** ss);
void   update_lda_seq_ss_serial(int time, const doc_t* doc, const lda_post* post, gsl_matrix** ss);
double fit_lda_seq_topics(lda_seq* model, gsl_matrix** ss, const int iter);
gsl_matrix* init_mu(const int ntopics, const int nmetafields);
void update_mu(const corpus_seq_t* data, lda_seq* model, const int iter); // Added by AG for rDIM

void write_lda_seq(const lda_seq* m, const char* root);
lda_seq* read_lda_seq(const char* root, corpus_seq_t* data);
void write_lda_seq_suffstats(lda_seq* m, gsl_matrix** topic_ss, const char* root);
lda_seq* new_lda_seq(corpus_seq_t* data, int W, int T, int K);
void free_lda_seq_model(lda_seq* m); // added by YH
void make_lda_from_seq_slice(lda* lda_m, lda_seq* lda_seq_m, int time);
static gsl_matrix* g_alloc(lda_seq* model, const corpus_seq_t* data, int time);

#ifdef SPARSE
// Some helper functions for spmatrix:
void gsl_spmatrix_get_row(gsl_vector* v, const gsl_spmatrix* M, const size_t x);
void gsl_spmatrix_set_col(gsl_spmatrix* M, const size_t x, const gsl_vector* v);
double g_round(double x, unsigned int digits);

// The Sparse Linear Solver functions // Added by AG
void make_term1( const corpus_seq_t* data, lda_seq* model);
void make_term1_csr3( const corpus_seq_t* data, lda_seq* model);
int  SparseSolve(const gsl_spmatrix* M, const gsl_vector* b, gsl_vector* x, const int s);
void add_noise_scaling(gsl_matrix *dest, const double epsilon);
int  MKLSolve_lapacke( gsl_matrix* gsl_A, gsl_matrix* gsl_b, gsl_matrix* gsl_x, 
		       mkl_factor* factor, bool& try_factor);
#endif

#endif
