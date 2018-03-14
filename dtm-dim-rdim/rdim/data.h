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

#ifndef DATA_H
#define DATA_H

#include "gsl-wrappers.h"
#include "param.h"
#include <stdio.h>
#include <stdlib.h>
#include <utility>

#ifdef __cplusplus 
extern "C" { 
#endif 

#ifdef SPARSE
#include <gsl/gsl_spmatrix.h>
#endif

#ifdef __cplusplus
} 
#endif

#ifdef __cplusplus 
extern "C" { 
#endif 
namespace mklspace {

#include "mkl_lapacke.h"
#include "mkl.h"

}
#ifdef __cplusplus 
} 
#endif

#define OFFSET 0

// Create the scaled beta distribution, which describes
// how much weight documents have after n years.
const int kScaledInfluenceMax = 200;

// This mean and variance are relative to the interval [0, 1.0].
const double kScaledInfluenceMean     =   10.0 / kScaledInfluenceMax;
const double kScaledInfluenceVariance = ((10.0 / kScaledInfluenceMax) * (10.0 / kScaledInfluenceMax));

/*
 * a document is a collection of counts and terms
 *
 */
typedef struct doc_t {
    int     total;
    int     nterms;
    int*    word;
    int*    count;
    double* lambda;        // A parameter for finding phi.
    double  log_likelihood; // Used for measuring perplexity.
    double* log_likelihoods;
} doc_t;

/*
 * a corpus is a collection of documents
 *
 */
typedef struct corpus_t {
    doc_t** doc;
    int     ndocs;
    int     nterms;
    int     max_unique;  // maximum number of unique terms in a document
} corpus_t;

/*
 * a sequence is a sequence of corpora
 *
 */
typedef struct corpus_seq_t {
#ifdef SPARSE
  gsl_spmatrix** tau; // |D_t| x nmetafiels for every t // AG TODO: This should really be in struct corpus
#else
  gsl_matrix**   tau; // |D_t| x nmetafiels for every t // AG TODO: This should really be in struct corpus
#endif
  corpus_t**     corpus;
  int            nterms;
  int            max_nterms;
  int            len;
  int            ndocs;
  int            nmetafields; // Added by AG
} corpus_seq_t;

typedef struct inf_var {
  gsl_matrix** doc_weights;              // T matrices of document weights. each matrix is |D_t| x K.
  gsl_matrix** renormalized_doc_weights; // T matrices of document weights. each matrix is |D_t| x K.
  int          ntime;
} inf_var;

/*
 * variational posterior structure
 *
 */
typedef struct sslm_var {
  // properties
  int W; // vocabulary size ~
  int T; // sequence length ~

  // variational parameters
  gsl_matrix* obs;             // observations, W x T ~
  
  // biproducts of the variational parameters
  double      obs_variance;    // observation variance
  double      chain_variance;  // chain variance
  gsl_vector* zeta;            // extra variational parameter, T
  gsl_matrix* e_log_prob;      // E log prob(w | t), W x T
  
  // convenient quantities for inference
  gsl_matrix* fwd_mean;       // forward posterior mean, W x T
  gsl_matrix* fwd_variance;   // forward posterior variance, W x T
  gsl_matrix* mean;           // posterior mean, W x T
  gsl_matrix* variance;       // posterior variance, W x T

  // These three were never used:
  //  gsl_matrix* mean_t;         // W x T
  //  gsl_matrix* variance_t;
  //  gsl_matrix* influence_sum_lgl;// The sum exp * w_phi_l
  
// Recent copy of w_phi_l.
  gsl_matrix* w_phi_l;          // W x T
  gsl_matrix* w_phi_sum;        // W x T
  gsl_matrix* w_phi_l_sq;       // Square term involving various
  gsl_matrix* m_update_coeff;   // Terms involving squares of W, l, and phi.
  gsl_matrix* m_update_coeff_g; // \sum_i=0..t phi_l(t) r(i-t)
  gsl_vector* T_vct;            // useful temporary vector
} sslm_var;

#ifdef MKL
typedef struct mkl_factor {
  MKL_INT*    ipiv;
  gsl_matrix* af;
  double*     r; // n-length array
  double*     c; // ibid
  int         n;
} mkl_factor;
#endif

typedef struct lda_seq {
  gsl_matrix*  eq32_term1; // |S| x |S| ~
  int          eq32_term1_nz; // Number of non-zeros is eq32_term1
  int          ntopics;    // number of topics ~ 
  int          nterms;     // number of terms ~
  int          nseq;       // length of sequence ~
  gsl_vector*  alpha;      // dirichlet parameters ~
  sslm_var**   topic;      // topic chains.
  inf_var*     influence;  // document weights
  gsl_matrix** influence_sum_lgl; // Sum of document weights at time t (see g in the regression formula)
  gsl_matrix*  mu;         // K x nmetafields // Added by AG 

  // These aren't used anywhere // AG
  //std::pair<int, float>**** top_doc_phis;  // T x D_t x n of document phis.

#ifdef MKL
  bool A_factored;
  bool try_factor;
  mkl_factor* A_factor; // TO CHECKPOINT
#endif
} lda_seq;

/*
 * functions
 *
 */
corpus_t*     read_corpus(           const char*         name);
corpus_seq_t* read_corpus_seq(       const char*         name, const int ntopics, const int nmetafields);
int           compute_max_nterms(    const corpus_seq_t* c);
gsl_matrix *  compute_total_counts(  const corpus_seq_t* c);
corpus_seq_t* make_seq_corpus_subset(corpus_seq_t*       all, int   start,    int end);
void          write_corpus(          corpus_t*           c,   char* filename);
void          write_corpus_seq(      corpus_seq_t*       c,   char* name);
corpus_seq_t* make_corpus_seq_subset(corpus_seq_t*       all, int   start,    int end);
corpus_t*     collapse_corpus_seq(   corpus_seq_t*       c);
double*       NewScaledInfluence(    int                 size);

#ifdef SPARSE
gsl_spmatrix* read_tau(              const char*         root,  const int t, 
				     const int           ndocs, const int nmetafields);
#else
gsl_matrix*   read_tau(              const char*         root,  const int t, 
				     const int           ndocs, const int nmetafields);
#endif

void free_corpus(corpus_t* data); // added by YH
void free_corpus_seq(corpus_seq_t* data_seq);  // added by YH
#endif
