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

#include "gflags.h"
#include "lda.h"
#include "omp.h"

int LDA_INFERENCE_MAX_ITER=25;

// Used for the phi update: how small can the smallest phi be?
// We assert around 1e-10, since that is around the smallest beta for a term.
static const double  kSmallLogNumber = -100.0;
static const double  kSmallNumber = exp(kSmallLogNumber);
static       double* scaled_influence = NULL;

DEFINE_double(lambda_convergence, 0.01, "Specifies the level of convergence required for lambda in the phi updates.");
DEFINE_bool(minimize_phi_iterations, false, "Limit number of iterations for phi-updates during dynamic fit to (CURRENT_EM_ITERATION+1)*2. The idea: A good estimate of lambda for the phi-updates is less important until the topic chains converge themselves.");
DECLARE_string(normalize_docs);
DECLARE_int32(max_number_time_points);
DECLARE_string(model);
DECLARE_string(mode);
DECLARE_double(sigma_d);
DECLARE_double(sigma_l);
DECLARE_int32(forward_window);
DECLARE_int32(threads);
DECLARE_int64(rng_seed);
DECLARE_int32(debug);
DECLARE_double(ntopics);
DECLARE_int32(kthreads);

/*
 * posterior inference for lda
 * time and doc_number are only necessary if
 * var is not NULL.
 */
double fit_lda_post(int         doc_number, int         time,      lda_post* p, lda_seq* var, gsl_matrix* g,
		    gsl_matrix* g3_matrix,  gsl_matrix* g4_matrix, gsl_matrix* g5_matrix) {
  init_lda_post(p);

  gsl_vector_view topic_view;
  gsl_vector_view renormalized_topic_view;

  if ((FLAGS_model == "fixed" || FLAGS_model == "rdim") && var && var->influence) {
    // Make sure this stays in scope while the posterior is in use!
    topic_view = gsl_matrix_row( var->influence->doc_weights[time], doc_number);
    renormalized_topic_view = gsl_matrix_row( var->influence->renormalized_doc_weights[time], doc_number);
    p->doc_weight = &topic_view.vector;
    p->renormalized_doc_weight = &renormalized_topic_view.vector;
  }

  double lhood     = compute_lda_lhood(p);
  double lhood_old = 0;
  double converged = 0;
  int iter = 0;

  do {
    iter++;
    lhood_old = lhood;
    update_gamma(p);

    if ((FLAGS_model == "fixed"  || FLAGS_model == "rdim") && var != NULL) {
      update_phi_fixed(doc_number, time, p, var, g3_matrix, g4_matrix, g5_matrix);
    } else if (FLAGS_model == "dtm" || var == NULL) {
      update_phi(doc_number, time, p, var, g);
    } else {
      outlog("Error. Unknonwn model '%s'.\n", FLAGS_model.c_str());
      exit(1);
    }

    lhood     = compute_lda_lhood(p);
    converged = fabs((lhood_old - lhood) / (lhood_old * p->doc->total));
  } while ((converged > LDA_INFERENCE_CONVERGED) && (iter <= LDA_INFERENCE_MAX_ITER));

  return(lhood);
}

// This assumes model == rdim or dim
double fit_lda_post_unrolled(int seqiter, int doc_number, int time, lda_post* p,
			     lda_seq* var, gsl_matrix* g, gsl_matrix* g3_matrix,
			     gsl_matrix* g4_matrix, gsl_matrix* g5_matrix) {
  const  int     K = p->model->ntopics, N = p->doc->nterms;
  const  double v2 = 1.0 / K;
  const  double v1 = ((double) p->doc->total) / K;
  const  double sigma_l2 = FLAGS_sigma_l * FLAGS_sigma_l;
  const  double sigma_d2 = FLAGS_sigma_d * FLAGS_sigma_d;
  const  double neg_sigma_l2      = (-1.0 * sigma_l2);
  int    max_internal_iter = LDA_INFERENCE_MAX_ITER+1;

  if (FLAGS_minimize_phi_iterations)
    max_internal_iter = (seqiter+1) * 2;

  gsl_matrix_set_all(p->phi,v2);

  for (int k = 0; k < K; k++)
    vset(p->gamma, k, vget(p->model->alpha,k) + v1);

  p->doc_weight = NULL;

  gsl_vector_view topic_view;
  gsl_vector_view renormalized_topic_view;
  
  if (var && var->influence) {
    topic_view = gsl_matrix_row( var->influence->doc_weights[time], doc_number);
    renormalized_topic_view = gsl_matrix_row( var->influence->renormalized_doc_weights[time], doc_number);
    p->doc_weight = &topic_view.vector;
    p->renormalized_doc_weight = &renormalized_topic_view.vector;
  }

  double gamma_sum = sum(p->gamma);
  double lhood = gsl_sf_lngamma(sum(p->model->alpha)) - gsl_sf_lngamma(gamma_sum);  
  vset(p->lhood, K, lhood);
  double digsum = gsl_sf_psi(gamma_sum);

  double influence_term = 0.0;

  for (int k = 0; k < K; k++) {
    if (p->doc_weight != NULL) {
      
      double influence_topic = gsl_vector_get(p->doc_weight, k);
      influence_term = - ((influence_topic * influence_topic + sigma_l2) * .5 / (sigma_d2));
    }
    
    double e_log_theta_k =  gsl_sf_psi(vget(p->gamma, k)) - digsum;
    double lhood_term    = (vget(p->model->alpha, k)-vget(p->gamma, k)) * e_log_theta_k +
                            gsl_sf_lngamma(vget(p->gamma, k)) - gsl_sf_lngamma(vget(p->model->alpha, k));
    
    for (int n = 0; n < N; n++) {
      lhood_term +=  p->doc->count[n] * /* v2 * */ (e_log_theta_k
	       	   + mget(p->model->topics, p->doc->word[n], k) - mget(p->log_phi, n, k));
    }
    lhood *= v2;
    
    vset(p->lhood, k, lhood_term);
    lhood = lhood + lhood_term + lhood + influence_term;
  }


  double lhood_old = 0;
  double converged = 0;
  int iter = 0;

  do {
    iter++;
    lhood_old = lhood;

    gsl_vector_memcpy(p->gamma, p->model->alpha);
    
    if (var == NULL) {
       double v;
       for (int k = 0; k < K; k++) {
         for (int n = 0; n < N; n++) {     
           v  = gsl_matrix_get(p->phi, n, k) * p->doc->count[n];
	   vinc(p->gamma, k, v);
         }
       }

       update_phi(doc_number, time, p, var, g);
    }
    else {
      if (scaled_influence == NULL) 
	scaled_influence = NewScaledInfluence(FLAGS_max_number_time_points);

      double dig[p->model->ntopics];
      double k_sum = 0.0;
      double v;

      for (int k = 0; k < K; k++) {
	for (int n = 0; n < N; n++) {     
	  v  = gsl_matrix_get(p->phi, n, k) * p->doc->count[n];
	  vinc(p->gamma, k, v);
	}
      }

      double dig_sum ;
      gsl_vector_view document_weights;
      double gamma_k;
	
      for (int k = 0; k < K; k++) {
        gamma_k = vget(p->gamma, k);
	dig[k] = gsl_sf_psi(gamma_k);
	k_sum  = k_sum + gamma_k;
      }

      dig_sum = gsl_sf_psi(k_sum);

      if (var && var->influence)
	document_weights = gsl_matrix_row(var->influence->doc_weights[time], doc_number);

#pragma omp parallel num_threads(FLAGS_kthreads)
{  
      if (time < var->nseq) {
#pragma omp for schedule(static)  
	for (int n=0; n < N; ++n) {
	  //if (n ==0) {
             //debug("OMP DEBUG A (%i threads)\n",omp_get_num_threads());
          //}

	  int w = p->doc->word[n];
	  double term_weight  = ((double) p->doc->count[n] / (double) p->doc->total);
	  double t2           = term_weight * term_weight;

	  for (int k = 0; k < K; ++k) {
	    double total,
	      doc_weight = gsl_vector_get(&document_weights.vector, k),
	      td = term_weight * doc_weight,
	      g3 = mget(g3_matrix, w, k),
	      g4 = mget(g4_matrix, w, k),
	      g5 = mget(g5_matrix, w, k);
	    
	    const double chain_variance = var->topic[k]->chain_variance;
	
	    total = (dig[k] + mget(p->model->topics, w, k))
	      - (g3 * td / chain_variance) - (td * g4 / chain_variance)
	      + (t2 * (-1.0 * (doc_weight * doc_weight + sigma_l2)) * g5 / chain_variance);
	    
	    mset(p->log_phi, n, k, total);
	  }
    
	  // Normalize in log space.
	  gsl_vector log_phi_row = gsl_matrix_row(p->log_phi, n).vector;
	  gsl_vector phi_row     = gsl_matrix_row(p->phi, n).vector;
	  log_normalize(&log_phi_row);
      
	  for (int i = 0; i < K; i++) 
	    vset(&phi_row, i, exp(vget(&log_phi_row, i)));
	}
      }
      else {
 #pragma omp for schedule(static)
	for (int n=0; n < N; ++n) {
	  int w = p->doc->word[n];
	  double term_weight  = ((double) p->doc->count[n] / (double) p->doc->total);
	  double t2           = term_weight * term_weight * neg_sigma_l2;
	  
	  for (int k = 0; k < K; ++k) {	
 	    double total,
	      g3 = mget(g3_matrix, w, k),
	      g4 = mget(g4_matrix, w, k),
	      g5 = mget(g5_matrix, w, k);
	    
	    const double chain_variance = var->topic[k]->chain_variance;	    
	    total = (dig[k] + mget(p->model->topics, w, k)) + (t2 * g5 / chain_variance);

	    mset(p->log_phi, n, k, total);
          }
   
	  // Normalize in log space.
	  gsl_vector log_phi_row = gsl_matrix_row(p->log_phi, n).vector;
	  gsl_vector phi_row     = gsl_matrix_row(p->phi, n).vector;
	  log_normalize(&log_phi_row);
       
	  for (int i = 0; i < K; i++) 
	    vset(&phi_row, i, exp(vget(&log_phi_row, i)));
        }
      }
    } // GOMP
    } // VAR IS NOT NULL _fixed_ run

    gamma_sum = sum(p->gamma);
    lhood = gsl_sf_lngamma(sum(p->model->alpha)) - gsl_sf_lngamma(gamma_sum);  
    vset(p->lhood, K, lhood);
    digsum = gsl_sf_psi(gamma_sum);
  
//#pragma omp parallel num_threads(FLAGS_kthreads)
    //{ 
    double influence_term = 0.0;
//#pragma omp for reduction(+:lhood) schedule(static) 
    for (int k = 0; k < K; k++) {
      if (p->doc_weight != NULL) {
      
	double influence_topic = gsl_vector_get(p->doc_weight, k);
	influence_term = - ((influence_topic * influence_topic + sigma_l2) * .5 / (sigma_d2));
      }
    
      double e_log_theta_k = gsl_sf_psi(vget(p->gamma, k)) - digsum;
      double lhood_term    = (vget(p->model->alpha, k)-vget(p->gamma, k)) * e_log_theta_k
	+ gsl_sf_lngamma(vget(p->gamma, k)) - gsl_sf_lngamma(vget(p->model->alpha, k));
    
      for (int n = 0; n < N; n++) {
	double p_phi = mget(p->phi, n, k);
	if (p_phi > 0)
	  lhood_term += p->doc->count[n] * p_phi * (e_log_theta_k
	  + mget(p->model->topics, p->doc->word[n], k) - mget(p->log_phi, n, k));
      }
    
      vset(p->lhood, k, lhood_term);

      lhood = lhood + lhood_term + lhood + influence_term;
    }
    //} // GOMP

    converged = fabs((lhood_old - lhood) / (lhood_old * p->doc->total));
  } while ((converged > LDA_INFERENCE_CONVERGED) &&
	   (iter <= LDA_INFERENCE_MAX_ITER)      &&
	    iter <  max_internal_iter);

  debug("PHI UPDATE: Used %d iterations on static LDA.\n", iter);
  
  return(lhood);
}

// This assume model == rdim or dim
double fit_lda_post_unrolled_serial_prenull(int doc_number, int time,
					    lda_post* p, lda_seq* var, gsl_matrix* g,
					    gsl_matrix* g3_matrix,  gsl_matrix* g4_matrix, 
					    gsl_matrix* g5_matrix) {
  const  int     K = p->model->ntopics, N = p->doc->nterms;
  const  double v2 = 1.0 / K;
  const  double v1 = ((double) p->doc->total) / K;
  const  double sigma_l2 = FLAGS_sigma_l * FLAGS_sigma_l;
  const  double sigma_d2 = FLAGS_sigma_d * FLAGS_sigma_d;
  const  double neg_sigma_l2 = (-1.0 * sigma_l2);
  
  gsl_matrix_set_all(p->phi,v2);

  for (int k = 0; k < K; k++)
    vset(p->gamma, k, vget(p->model->alpha,k) + v1);

  p->doc_weight = NULL;

  gsl_vector_view topic_view;
  gsl_vector_view renormalized_topic_view;

  double gamma_sum = sum(p->gamma);
  double lhood = gsl_sf_lngamma(sum(p->model->alpha)) - gsl_sf_lngamma(gamma_sum);  
  vset(p->lhood, K, lhood);
  double digsum = gsl_sf_psi(gamma_sum);
  double influence_term = 0.0;
  
  for (int k = 0; k < K; k++) {
    if (p->doc_weight != NULL) {
      double influence_topic = gsl_vector_get(p->doc_weight, k);
      influence_term = - ((influence_topic * influence_topic + sigma_l2) * .5 / (sigma_d2));
    }
    
    double e_log_theta_k = gsl_sf_psi(vget(p->gamma, k)) - digsum;
    double lhood_term    = (vget(p->model->alpha, k)-vget(p->gamma, k)) * e_log_theta_k +
                            gsl_sf_lngamma(vget(p->gamma, k)) - gsl_sf_lngamma(vget(p->model->alpha, k));
    
    for (int n = 0; n < N; n++) {
      lhood_term += p->doc->count[n] * (e_log_theta_k
	       	  + mget(p->model->topics, p->doc->word[n], k) 
		  - mget(p->log_phi, n, k));
    }

    lhood *= v2;
    
    vset(p->lhood, k, lhood_term);
    lhood = lhood + lhood_term + lhood + influence_term;
  }

  double lhood_old = 0;
  double converged = 0;
  int iter = 0;

  do {
    iter++;
    lhood_old = lhood;

    gsl_vector_memcpy(p->gamma, p->model->alpha);
    
    double v;
    for (int k = 0; k < K; k++) {
      for (int n = 0; n < N; n++) {     
	v  = gsl_matrix_get(p->phi, n, k) * p->doc->count[n];
	vinc(p->gamma, k, v);
      }
    }

    update_phi(doc_number, time, p, var, g);

    gamma_sum = sum(p->gamma);
    lhood = gsl_sf_lngamma(sum(p->model->alpha)) - gsl_sf_lngamma(gamma_sum);  
    vset(p->lhood, K, lhood);
    digsum = gsl_sf_psi(gamma_sum);
  
    double influence_term = 0.0;
    for (int k = 0; k < K; k++) {
      if (p->doc_weight != NULL) {
      
	double influence_topic = gsl_vector_get(p->doc_weight, k);
	influence_term = - ((influence_topic * influence_topic + sigma_l2) * .5 / (sigma_d2));
      }
    
      double e_log_theta_k = gsl_sf_psi(vget(p->gamma, k)) - digsum;
      double lhood_term    = (vget(p->model->alpha, k)-vget(p->gamma, k)) * e_log_theta_k
	+ gsl_sf_lngamma(vget(p->gamma, k)) - gsl_sf_lngamma(vget(p->model->alpha, k));
    
      for (int n = 0; n < N; n++) {
	double p_phi = mget(p->phi, n, k);
	
	if (p_phi > 0)
	  lhood_term += p->doc->count[n] * p_phi * (e_log_theta_k
	              + mget(p->model->topics, p->doc->word[n], k) - mget(p->log_phi, n, k));
      }
      
      vset(p->lhood, k, lhood_term);
      lhood = lhood + lhood_term + lhood + influence_term;
    }
    
    converged = fabs((lhood_old - lhood) / (lhood_old * p->doc->total));
  } while ((converged > LDA_INFERENCE_CONVERGED) && (iter <= LDA_INFERENCE_MAX_ITER));

  return(lhood);
}

/*
 * initialize variational posterior
 *
 */
void init_lda_post(lda_post* p) {
  int K = p->model->ntopics, N = p->doc->nterms;
  double v2 = 1.0/K;
  double v1 = ((double) p->doc->total) / K;

  gsl_matrix_set_all(p->phi, v2);

  for (int k = 0; k < K; k++)
    vset(p->gamma, k, vget(p->model->alpha,k) + v1);

  p->doc_weight = NULL;
}

/*
 * update variational dirichlet parameters
 *
 */
void update_gamma(lda_post* p) {
  int K = p->model->ntopics, N = p->doc->nterms;
  gsl_vector_memcpy(p->gamma, p->model->alpha);

#pragma omp parallel num_threads(FLAGS_kthreads)
{
  double v;
#pragma omp for schedule(static) 
  for (int k = 0; k < K; k++) {
    if (k ==0) {
      debug("OMP DEBUG C (%i threads)\n",omp_get_num_threads());
    }
    for (int n = 0; n < N; n++) {     
      v  = gsl_matrix_get(p->phi, n, k);
      v *= p->doc->count[n];
      vinc(p->gamma, k, v);
    }
  }
  }//GOMP
}

/*
 * update variational multinomial parameters
 *
 */
void update_phi(int doc_number, int time, lda_post* p, lda_seq* var, gsl_matrix* g) {
  int i, k, n, K = p->model->ntopics, N = p->doc->nterms;
  double dig[p->model->ntopics];

  for (k = 0; k < K; k++) 
    dig[k] = gsl_sf_psi(vget(p->gamma, k));

  for (n = 0; n < N; n++) {
    int w = p->doc->word[n];

    for (k = 0; k < K; k++) 
      mset(p->log_phi, n, k, dig[k] + mget(p->model->topics, w, k));

    gsl_vector log_phi_row = gsl_matrix_row(p->log_phi, n).vector;
    gsl_vector phi_row = gsl_matrix_row(p->phi, n).vector;
    log_normalize(&log_phi_row);

    for (i = 0; i < K; i++) 
      vset(&phi_row, i, exp(vget(&log_phi_row, i)));
  }
}

void update_phi_fixed(int doc_number, int time, lda_post*   p,         lda_seq*    var,
		      gsl_matrix*    g3_matrix, gsl_matrix* g4_matrix, gsl_matrix* g5_matrix) {
  // Hate to do this, but I had problems allocating this data structure. // SG
  if (scaled_influence == NULL) 
    scaled_influence = NewScaledInfluence(FLAGS_max_number_time_points);

  int 
    K = p->model->ntopics, 
    N = p->doc->nterms;
  double dig[p->model->ntopics];
  double k_sum = 0.0;
  double sigma_l2 = FLAGS_sigma_l * FLAGS_sigma_l;

//#pragma omp parallel default(shared) num_threads(FLAGS_kthreads) // Old way
#pragma omp parallel default(shared) num_threads(FLAGS_threads)
  {
  double dig_sum ;
  gsl_vector_view document_weights;
  double gamma_k;

  #pragma omp for schedule(static) reduction(+:k_sum)
  for (int k = 0; k < K; k++) {
    if (k ==0) {
      debug("OMP DEBUG C (%i threads)\n",omp_get_num_threads());
    }
    gamma_k = vget(p->gamma, k);
    dig[k] = gsl_sf_psi(gamma_k);
    k_sum  = k_sum + gamma_k;
  }

  dig_sum = gsl_sf_psi(k_sum);

  if (var && var->influence)
    document_weights = gsl_matrix_row(var->influence->doc_weights[time], doc_number);

  if (time < var->nseq) {
  #pragma omp for schedule(static)
    for (int n=0; n < N; ++n) {
      if (n ==0) {
        debug("OMP DEBUG D (%i threads)\n",omp_get_num_threads());
      }
      int w = p->doc->word[n];
      double term_weight  = ((double) p->doc->count[n] / (double) p->doc->total);
      double t2           = term_weight * term_weight;

      for (int k = 0; k < K; ++k) {
	double total,
	  doc_weight = gsl_vector_get(&document_weights.vector, k),
	  td = term_weight * doc_weight,
	  g3 = mget(g3_matrix, w, k),
	  g4 = mget(g4_matrix, w, k),
	  g5 = mget(g5_matrix, w, k);
	
	const double chain_variance = var->topic[k]->chain_variance;
	
	total = (dig[k] + mget(p->model->topics, w, k))
	  - (g3 * td / chain_variance) - (td * g4 / chain_variance)
	  + (t2 * (-1.0 * (doc_weight * doc_weight + sigma_l2)) * g5 / chain_variance);

	mset(p->log_phi, n, k, total);
      }
    
      // Normalize in log space.
      gsl_vector log_phi_row = gsl_matrix_row(p->log_phi, n).vector;
      gsl_vector phi_row     = gsl_matrix_row(p->phi, n).vector;
      log_normalize(&log_phi_row);
      
      for (int i = 0; i < K; i++) 
	vset(&phi_row, i, exp(vget(&log_phi_row, i)));
    }
  }
  else {
  #pragma omp for schedule(static)
    for (int n=0; n < N; ++n) {
	if (n ==0) {
          debug("OMP DEBUG E (%i threads)\n",omp_get_num_threads());
        }
	int w = p->doc->word[n];
	double term_weight  = ((double) p->doc->count[n] / (double) p->doc->total);
	double t2           = term_weight * term_weight;

	for (int k = 0; k < K; ++k) {	
	  double total,
	    g3 = mget(g3_matrix, w, k),
	    g4 = mget(g4_matrix, w, k),
	    g5 = mget(g5_matrix, w, k);

	  const double chain_variance = var->topic[k]->chain_variance;
       
	  total = (dig[k] + mget(p->model->topics, w, k))
	    + (t2 * (-1.0 * sigma_l2) * g5 / chain_variance); 

	  mset(p->log_phi, n, k, total);
	}
   
	// Normalize in log space.
	gsl_vector log_phi_row = gsl_matrix_row(p->log_phi, n).vector;
	gsl_vector phi_row     = gsl_matrix_row(p->phi, n).vector;
	log_normalize(&log_phi_row);
      
	for (int i = 0; i < K; i++) 
	  vset(&phi_row, i, exp(vget(&log_phi_row, i)));
    }
  }
  } // GOMP
}

/*
 * comput the likelihood bound
 */
double compute_lda_lhood(lda_post* p) {
  int K = p->model->ntopics, N = p->doc->nterms;

  double gamma_sum = sum(p->gamma);
  double lhood = gsl_sf_lngamma(sum(p->model->alpha)) - gsl_sf_lngamma(gamma_sum);
  
  vset(p->lhood, K, lhood);
  double digsum = gsl_sf_psi(gamma_sum);
  double sigma_l2 = FLAGS_sigma_l * FLAGS_sigma_l;
  double sigma_d2 = FLAGS_sigma_d * FLAGS_sigma_d;
  bool influence_model = (FLAGS_model == "rdim" || FLAGS_model == "fixed" || FLAGS_model == "dim") ;
  
#pragma omp parallel num_threads(FLAGS_kthreads)
  {
  double influence_term = 0.0;
#pragma omp for reduction(+:lhood) schedule(static)
  for (int k = 0; k < K; k++) {
    if (K ==0) {
      debug("OMP DEBUG F (%i threads)\n",omp_get_num_threads());
    }
    if (p->doc_weight != NULL) {
      
      double influence_topic = gsl_vector_get(p->doc_weight, k);
      if (influence_model) 
	influence_term = - ((influence_topic * influence_topic + sigma_l2) * .5 / (sigma_d2));
    }
    
    double e_log_theta_k = gsl_sf_psi(vget(p->gamma, k)) - digsum;
    double lhood_term    = (vget(p->model->alpha, k)-vget(p->gamma, k)) * e_log_theta_k +
                            gsl_sf_lngamma(vget(p->gamma, k)) - gsl_sf_lngamma(vget(p->model->alpha, k));
    
    for (int n = 0; n < N; n++) {
      double p_phi = mget(p->phi, n, k);
      if (p_phi > 0)
	lhood_term += p->doc->count[n] * p_phi * (e_log_theta_k
	            + mget(p->model->topics, p->doc->word[n], k) - mget(p->log_phi, n, k));
    }
    
    vset(p->lhood, k, lhood_term);
    lhood = lhood + lhood_term + lhood + influence_term;
  }
  } //GOMP 

  return(lhood);
}

/*
 * compute expected sufficient statistics for a corpus
 *
 */
double lda_e_step(lda* model, corpus_t* data, lda_suff_stats* ss) {
  omp_set_dynamic(0);

  double lhood = 0;
  double etime = omp_get_wtime();
  
#pragma omp parallel shared(lhood) num_threads(FLAGS_threads)
  { 
    outlog("(LDA-init) we have %i GOMP threads (%i specified on CLI).\n",omp_get_num_threads(),FLAGS_threads);
    if (ss != NULL) 
      gsl_matrix_set_all(ss->topics_ss, 0.0);
    
    lda_post* p = new_lda_post(model->ntopics, data->max_unique);
    p->model    = model;

    for (int k = 0; k < model->ntopics; k++) {
      outlog("(LDA-init) [k=%i, gomp_thread=%i (%i)]\n", k,omp_get_thread_num(), omp_get_num_threads());
#pragma omp for reduction(+:lhood) schedule(static)
      for (int d = 0; d < data->ndocs; d++) {
	if (d == 0) {
           debug("OMP DEBUG H (%i threads)\n",omp_get_num_threads());
        }
	p->doc = data->doc[d];

	//lhood += fit_lda_post(d, 0, p, NULL, NULL, NULL, NULL, NULL);
	lhood = lhood + fit_lda_post_unrolled_serial_prenull(d, 0, p, NULL, NULL, NULL, NULL, NULL);
	
	//	if (ss != NULL) // This can never be null here:
	for (int n = 0; n < p->doc->nterms; n++)
	  minc(ss->topics_ss, p->doc->word[n], k, mget(p->phi, n, k) * p->doc->count[n]);
      } // d
    } // k

    free_lda_post(p); // added by YH
  } // End OMP 
  //omp_set_dynamic(1); // TODO: old way was dyn here

  debug("LDA INIT E-STEP took %f\n", omp_get_wtime()-etime);

  return(lhood);
}

/*
 * compute MLE topics from sufficient statistics
 *
 */
double lda_m_step(lda* model, lda_suff_stats* ss) {
  double lhood = 0;

  for (int k = 0; k < model->ntopics; k++) {
    outlog("(LDA) [k=%i] M-Step\n", k);

    gsl_vector ss_k  = gsl_matrix_column(ss->topics_ss, k).vector;
    gsl_vector log_p = gsl_matrix_column(model->topics, k).vector;

    if (LDA_USE_VAR_BAYES == 0) {
      gsl_blas_dcopy(&ss_k, &log_p);
      normalize(&log_p);
      vct_log(&log_p);
    }
    else {
      double digsum    = sum(&ss_k)+model->nterms*LDA_TOPIC_DIR_PARAM;
      digsum           = gsl_sf_psi(digsum);
      double param_sum = 0;
      
      for (int w = 0; w < model->nterms; w++) {
	double param    = vget(&ss_k, w) + LDA_TOPIC_DIR_PARAM;
	param_sum      += param;
	double elogprob = gsl_sf_psi(param) - digsum;
	vset(&log_p, w, elogprob);
	lhood          += (LDA_TOPIC_DIR_PARAM - param) * elogprob + gsl_sf_lngamma(param);
      }
      
      lhood -= gsl_sf_lngamma(param_sum);
    }
  }
  
  return(lhood);
}

/*
 * read sufficient statistics
 *
 */
void write_lda_suff_stats(lda_suff_stats* ss, char* name) {
  mtx_fprintf(name, ss->topics_ss);
}

lda_suff_stats* read_lda_suff_stats(char* filename, int ntopics, int nterms) {
  lda_suff_stats* ss = (lda_suff_stats*) malloc(sizeof(lda_suff_stats));
  ss->topics_ss = gsl_matrix_alloc(nterms, ntopics);
  mtx_fscanf(filename, ss->topics_ss);

  return(ss);
}

/*
 * new lda model and sufficient statistics
 *
 */
lda* new_lda_model(int ntopics, int nterms) {
  lda* m = (lda*) malloc(sizeof(lda));

  m->ntopics = ntopics;
  m->nterms  = nterms;
  m->topics  = gsl_matrix_calloc(nterms, ntopics);
  m->alpha   = gsl_vector_calloc(ntopics);

  return(m);
}

void free_lda_model(lda* m) {
  gsl_matrix_free(m->topics);
  gsl_vector_free(m->alpha);
  free(m);
}

lda_suff_stats* new_lda_suff_stats(lda* model) {
  lda_suff_stats* ss = (lda_suff_stats*) malloc(sizeof(lda_suff_stats));
  ss->topics_ss = gsl_matrix_calloc(model->nterms, model->ntopics);
  
  return(ss);
}

void reset_lda_suff_stats(lda_suff_stats* ss) {
  gsl_matrix_set_all(ss->topics_ss, 0.0);
}

lda_post* new_lda_post(int ntopics, int max_length) {
  lda_post* p = (lda_post*) malloc(sizeof(lda_post));
  p->phi      = gsl_matrix_calloc(max_length, ntopics);
  p->log_phi  = gsl_matrix_calloc(max_length, ntopics);
  p->gamma    = gsl_vector_calloc(ntopics);
  p->lhood    = gsl_vector_calloc(ntopics + 1);

  return(p);
}

void free_lda_post(lda_post* p) {
  gsl_matrix_free(p->phi);
  gsl_matrix_free(p->log_phi);
  gsl_vector_free(p->gamma);
  gsl_vector_free(p->lhood);
  free(p);
}

/*
 * initalize LDA SS from random
 *
 */
void initialize_lda_ss_from_random(corpus_t* data, lda_suff_stats* ss) {
  int k, n;
  gsl_rng * r = new_random_number_generator((long) FLAGS_rng_seed);

  for (k = 0; k < ss->topics_ss->size2; k++) {
    gsl_vector topic = gsl_matrix_column(ss->topics_ss, k).vector;
    gsl_vector_set_all(&topic, 0);

    for (n = 0; n < topic.size; n++) 
      vset(&topic, n, gsl_rng_uniform(r) + 0.5 / data->ndocs + 4.0);
  }

  gsl_rng_free(r); // added by YH
}

/*
 * initialize sufficient statistics from a document collection
 *
 */
void initialize_lda_ss_from_data(corpus_t* data, lda_suff_stats* ss) {
  int k, n, i, w;
  gsl_rng * r = new_random_number_generator((long) FLAGS_rng_seed);

  for (k = 0; k < ss->topics_ss->size2; k++) {
    gsl_vector topic = gsl_matrix_column(ss->topics_ss, k).vector;

    for (n = 0; n < LDA_SEED_INIT; n++) {
      int d = floor(gsl_rng_uniform(r) * data->ndocs);
      doc_t* doc = data->doc[d];

      for (i = 0; i < doc->nterms; i++) 
	vinc(&topic, doc->word[n], doc->count[n]);
    }

    for (w = 0; w < topic.size; w++) 
      vinc(&topic, w, LDA_INIT_SMOOTH + gsl_rng_uniform(r));
  }

  gsl_rng_free(r); // Added by AG
}

/*
 * write LDA model
 *
 */
void write_lda(lda* model, char* name) {
  char filename[400];

  sprintf(filename, "%s.beta", name);
  mtx_fprintf(filename, model->topics);
  sprintf(filename, "%s.alpha", name);
  vct_fprintf(filename, model->alpha);
}

/*
 * read LDA
 *
 */
lda* read_lda(int ntopics, int nterms, char* name) {
  char filename[400];

  lda* model = new_lda_model(ntopics, nterms);
  sprintf(filename, "%s.beta", name);
  mtx_fscanf(filename, model->topics);
  sprintf(filename, "%s.alpha", name);
  vct_fscanf(filename, model->alpha);

  return(model);
}

void lda_em(lda* model, lda_suff_stats* ss, corpus_t* data, int max_iter, char* outname) {
  int    iter      = 0;
  double lhood     = lda_e_step(model, data, ss);
  double old_lhood = 0;
  double converged = 0;
  double m_lhood   = lda_m_step(model, ss);
  
  outlog( "(LDA) Initial likelihood = %f\n", lhood);

  do {
    iter++;
    old_lhood      = lhood;
    double e_lhood = lda_e_step(model, data, ss);
	
    m_lhood   =  lda_m_step(model, ss);
    lhood     =  e_lhood + m_lhood;
    converged = (old_lhood - lhood) / (old_lhood);
     
    outlog("(LDA) Iteration   = %d (%d)\n", iter, max_iter);
    outlog("(LDA) Likelihood  = % 10.3f\n", lhood);
    outlog("(LDA) M, E lhood  = % 10.3f, % 10.3f\n", m_lhood, e_lhood);
    outlog("(LDA) convergence = % 5.3e\n", converged);
  } while (((converged > LDA_EM_CONVERGED) || (iter <= 5)) && (iter < max_iter));
  
  write_lda(model, outname);
  return;
}
