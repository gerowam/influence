// Authors: Aaron Gerow (gerow@uchicago.edu)
//
// Copyright 2016 Aaron Gerow
// All Rights Reserved.
//
// rDIM version: Copyright 2015 Aaron Gerow, Yuening Hu, Jordan Boyd-Graber, James Evans and David Blei
// All Rights Reserved.
//
// See the README for this package for details about modifying or
// distributing this software.

#include <unistd.h>
#include "checkpoint.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "gflags.h"

DECLARE_int32(debug);

gsl_matrix* read_double_matrix(const char* descriptor, const char* dir) {
  char fname[BUFSIZ];
  gsl_matrix*   ret;
  hid_t         fid;
  herr_t        status;
  unsigned long s1, s2;
  
  sprintf(fname, "%s%s.h5", dir, descriptor);
  fid    = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
  H5LTget_attribute_ulong(fid, descriptor, "size1", &s1);
  H5LTget_attribute_ulong(fid, descriptor, "size2", &s2);

  if (s1 < 1 || s2 < 1) {
    outlog("ERROR: Got matrix->size1 < 1 for %s/$s (%i x %i), exiting...\n", dir, descriptor,s1,s2);
    _exit(24);
  }

  debug("Allocating %s to be read (%i x %i) ... ", fname, s1, s2);
  ret = gsl_matrix_alloc(s1, s2);

  status = H5LTread_dataset_double(fid, descriptor, ret->data);
  status = H5Fclose(fid);

  if (status == -1) {
    outlog("HDF5 returned -1 on reading %s/$s: %i\n", dir, descriptor, (int) status);
    _exit(22);
  }

  debug("read (%i x %i)\n", ret->size1, ret->size2);

  return ret;
}

gsl_vector* read_double_vector(const char* descriptor, const char* dir) {
  char fname[BUFSIZ];
  gsl_vector* ret;
  hid_t       fid;
  herr_t      status;
  int         s;
  
  sprintf(fname, "%s%s.h5", dir, descriptor);
  fid    = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
  H5LTget_attribute_int(  fid, descriptor, "size", &s);

  if (s < 1) {
    outlog("ERROR: Got vector->size < 1 for %s/$s (%i), exiting...\n", dir, descriptor,s);
    _exit(23);
  }
    
  ret    = gsl_vector_alloc(s);
  H5LTread_dataset_double(fid, descriptor, ret->data);
  status = H5Fclose(fid);

  if (status == -1) {
    outlog("HDF5 returned -1 on reading %s/$s: %i\n", dir, descriptor, (int) status);
    _exit(22);
  }

  return ret;
}

herr_t attach_double_matrix(const char* descriptor, const char* dir, gsl_matrix* mat) {
  if (mat != NULL && mat->data != NULL) {
    hid_t   fid;
    herr_t  status;
    hsize_t dims[2] = {mat->size1, mat->size2};
    char    fname[BUFSIZ];
    unsigned long 
      s1 = mat->size1,
      s2 = mat->size2;

    debug("CHECKP(attach): %s, (%i x %i)\n", descriptor, s1, s2);
    sprintf(fname, "%s%s.h5", dir, descriptor);
    fid    = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    status = H5LTmake_dataset(fid, descriptor, 2, dims, H5T_NATIVE_DOUBLE, mat->data);
    H5LTset_attribute_ulong(fid, descriptor, "size1", &s1, 1);
    H5LTset_attribute_ulong(fid, descriptor, "size2", &s2, 1);
  
    status = H5Fclose(fid);
  
    return status;
  }
  else
    debug("CHECKP(attach): %s points to NULL, skipping...\n", descriptor);
    return 1;
}

herr_t attach_double_vector(const char* descriptor, const char* dir, gsl_vector* vec) {
  if (vec != NULL && vec->data != NULL) {
    hid_t   fid;
    herr_t  status;
    hsize_t dims[1] = {vec->size};
    char    fname[BUFSIZ];
    unsigned long s = (int) vec->size;
    
    sprintf(fname, "%s%s.h5",dir, descriptor);
    fid    = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    status = H5LTmake_dataset(fid, descriptor, 1, dims, H5T_NATIVE_DOUBLE, vec->data);
    H5LTset_attribute_ulong(fid, descriptor, "size", &s, 1);
    
    status = H5Fclose(fid);
    
    return status;
  }
  else
    return 1;
}

// Needs _cplusplus
// breaks DTM readme and legacy Makefile:
int read_checkpoint(const char*   indir,            lda_seq    *&model,
		    gsl_matrix**  topic_suffstats,  gsl_matrix *&gammas,
		    gsl_matrix   *&lhoods,          gsl_matrix *&heldout_gammas,
		    gsl_matrix   *&heldout_lhoods,  double&     bound,
		    double&       heldout_bound,    double&     convergence,
		    double&       old_bound,        int&        iter,
		    short&        final_iters_flag, const bool  holdout,
 	            unsigned int& last_iter,        const int   K) { 

  herr_t status = 0;
  char   fname[BUFSIZ];
  char   rootk[BUFSIZ];
  char   roott[BUFSIZ];
  char   root[BUFSIZ];
  char   descriptor[BUFSIZ];
  hid_t  fid;
  int    attr_W, attr_T;
  double attr_obs_variance, attr_chain_variance;
  int    val;

  if (!directory_exist(indir)) {
    outlog("ERROR: Recovery director '%s' does not exist.", (char*) indir);
    return 2;
  }

  sprintf(root, "%s/", indir);

  debug("CHECKR 1%s\n","");
  sprintf(fname, "%s/attributes.h5", root); 
  fid    = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);

  H5LTget_attribute_int(   fid, "attrs", "model_ntopics",        &model->ntopics);
  H5LTget_attribute_int(   fid, "attrs", "model_nterms",         &model->nterms);
  H5LTget_attribute_int(   fid, "attrs", "model_nseq",           &model->nseq);
  H5LTget_attribute_int(   fid, "attrs", "model_influence_ntime",&model->influence->ntime);
  H5LTget_attribute_double(fid, "attrs", "bound",                &bound);
  H5LTget_attribute_double(fid, "attrs", "convergence",          &convergence);
  H5LTget_attribute_double(fid, "attrs", "old_bound",            &old_bound);
  H5LTget_attribute_int(   fid, "attrs", "iter",                 &iter);
  H5LTget_attribute_short( fid, "attrs", "final_iters_flag",     &final_iters_flag);
  H5LTget_attribute_uint(  fid, "attrs", "last_iter",            &last_iter);
  H5LTget_attribute_int(   fid, "attrs", "model_eq32_term1_nz", &model->eq32_term1_nz);
  H5LTget_attribute_int(   fid, "attrs", "model_A_factored",     &val);
  model->A_factored = val == 1 ? true : false;
  H5LTget_attribute_int(   fid, "attrs", "model_try_factor",     &val);
  model->try_factor = val == 1 ? true : false;

  debug("CHECKR 2a%s\n","");
  gammas = read_double_matrix("gammas", root);
  lhoods = read_double_matrix("lhoods", root);

  if (holdout) {
    debug("CHECKR 3%s\n","");
    H5LTget_attribute_double(fid, "attrs", "heldout_bound", &heldout_bound);
    heldout_gammas = read_double_matrix("heldout_gammas", root);
    heldout_lhoods = read_double_matrix("heldout_lhoods", root);
  }

  if (model->A_factored) {
    debug("CHECKP 4a: %i\n",(int) status);
    H5LTget_attribute_int(fid,"attrs","model_A_factor_n", &model->A_factor->n);
    model->A_factor->af = read_double_matrix("model_A_factor_af", root);

    gsl_vector *r    = read_double_vector("model_A_factor_r",    root);
    gsl_vector *c    = read_double_vector("model_A_factor_c",    root);
    gsl_vector *ipiv = read_double_vector("model_A_factor_ipiv", root);

    model->A_factor->ipiv  = (MKL_INT*) malloc(sizeof(MKL_INT) * model->A_factor->n);
    model->A_factor->r     = (double*)  malloc(sizeof(double)  * model->A_factor->n);
    model->A_factor->c     = (double*)  malloc(sizeof(double)  * model->A_factor->n);

    for (int i=0; i < model->A_factor->n; i++) { 
      model->A_factor->ipiv[i] = gsl_vector_get(ipiv, i);
      model->A_factor->r[i]    = gsl_vector_get(r,    i);
      model->A_factor->c[i]    = gsl_vector_get(c,    i);
    }
    
    gsl_vector_free(r);
    gsl_vector_free(c);
    gsl_vector_free(ipiv);
  }

  debug("CHECKR 4%s\n","");
  model->mu    = read_double_matrix("model_mu",    root);
  model->alpha = read_double_vector("model_alpha", root);

  for (int t=0; t < model->nseq; t++) {
    debug("CHECKR T%i\n", (int) t);    
    sprintf(roott, "%s/t%i-", indir, t);

    model->influence->doc_weights[t]              = read_double_matrix("model_influence_doc_weights", roott);
    model->influence->renormalized_doc_weights[t] = read_double_matrix("model_influence_renormalized_doc_weights", roott);
  }

  for (int k=0; k<K; k++) {
    debug("CHECKR K%i\n",  (int) k);

    sprintf(rootk, "%s/k%i-", root, (int) k);

    sprintf(descriptor, "model_topic_W-k%i", (int) k);
    H5LTget_attribute_int(fid, "attrs", descriptor, &attr_W);
    sprintf(descriptor, "model_topic_T-k%i", (int) k);
    H5LTget_attribute_int(fid, "attrs",    descriptor, &attr_T);
    sprintf(descriptor, "model_topic_obs_variance-k%i", (int) k);
    H5LTget_attribute_double(fid, "attrs", descriptor, &attr_obs_variance);
    sprintf(descriptor, "model_topic_chain_variance-k%i", (int) k);
    H5LTget_attribute_double(fid, "attrs", descriptor, &attr_chain_variance);

    if (attr_W != 0)
       model->topic[k]->W = attr_W;
    if (attr_T != 0)
       model->topic[k]->T = attr_T;
    if (attr_obs_variance != 0.0)
      model->topic[k]->obs_variance = attr_obs_variance;
    if (attr_chain_variance != 0.0)
      model->topic[k]->chain_variance = attr_chain_variance;

    model->topic[k]->obs               = read_double_matrix("model_topic_obs",               rootk);
    model->topic[k]->e_log_prob        = read_double_matrix("model_topic_e_log_prob",        rootk);
    model->topic[k]->mean              = read_double_matrix("model_topic_mean",              rootk);
    model->topic[k]->variance          = read_double_matrix("model_topic_variance",          rootk);
    // The ones that get free()'d but not NULLed at each M step:
    //    model->topic[k]->mean_t            = read_double_matrix("model_topic_mean_t",            rootk);
    //    model->topic[k]->variance_t        = read_double_matrix("model_topic_variance_t",        rootk);
    model->topic[k]->w_phi_l           = read_double_matrix("model_topic_w_phi_l",           rootk);
    model->topic[k]->w_phi_sum         = read_double_matrix("model_topic_w_phi_sum",         rootk);
    model->topic[k]->w_phi_l_sq        = read_double_matrix("model_topic_w_phi_l_sq",        rootk);
    model->topic[k]->m_update_coeff    = read_double_matrix("model_topic_m_update_coeff",   rootk);
    model->topic[k]->m_update_coeff_g  = read_double_matrix("model_topic_m_update_coeff_g", rootk);
    model->influence_sum_lgl[k]        = read_double_matrix("model_topic_influence_sum_lgl", rootk);
    topic_suffstats[k]                 = read_double_matrix("topic_suffstats",               rootk);
  }
  
  debug("CHECKR 8\n","");
  status = H5Fclose(fid);

  return (int) status;
}

int write_checkpoint(const char*  outdir,           lda_seq*     model,
		     gsl_matrix** topic_suffstats,  gsl_matrix*  gammas,
		     gsl_matrix*  lhoods,           gsl_matrix*  heldout_gammas,
		     gsl_matrix*  heldout_lhoods,   const double bound,
		     const double heldout_bound,    const double convergence,
		     const double old_bound,        const int    iter,
		     const short  final_iters_flag, const bool   holdout,
		     const unsigned int last_iter,  const int    K) {

  hid_t   fid,did,sid;
  herr_t  status;
  char    fname[BUFSIZ];
  char    root[BUFSIZ];
  char    rootk[BUFSIZ];
  char    roott[BUFSIZ];
  char    descriptor[BUFSIZ];
  hsize_t dims[1] = { 12 + (K * 4) };
  int     val;

  // Make the directory, if it doesn't exist:
  if (!directory_exist(outdir))
    make_directory((char*) outdir);

  sprintf(root, "%s/iter-%i", outdir, iter);
  if (!directory_exist(root))
    make_directory((char*) root);

  sprintf(root, "%s/", root);
  
  debug("CHECKP 1: %i\n",(int) status);
  // Attributes:
  sprintf(fname, "%s/attributes.h5",root, iter);
  fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  sid = H5Screate_simple(1, dims, NULL);
  did = H5Dcreate2(fid, "attrs", H5T_NATIVE_INT, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dclose(did);
  H5Sclose(sid);

  debug("CHECKP 2: %i\n",(int) status);
  
  // Model attributes attributed
  H5LTset_attribute_int(   fid,"attrs","model_ntopics",        &model->ntopics,          1);
  H5LTset_attribute_int(   fid,"attrs","model_nterms",         &model->nterms,           1);
  H5LTset_attribute_int(   fid,"attrs","model_nseq",           &model->nseq,             1);
  H5LTset_attribute_int(   fid,"attrs","model_influence_ntime",&model->influence->ntime, 1);
  H5LTset_attribute_double(fid,"attrs","bound",                &bound,                   1);
  H5LTset_attribute_double(fid,"attrs","convergence",          &convergence,             1);
  H5LTset_attribute_double(fid,"attrs","old_bound",            &old_bound,               1);
  H5LTset_attribute_int(   fid,"attrs","iter",                 &iter,                    1);
  H5LTset_attribute_short( fid,"attrs","final_iters_flag",     &final_iters_flag,        1);
  H5LTset_attribute_uint(  fid,"attrs","last_iter",            &last_iter,               1);
  H5LTset_attribute_int(   fid,"attrs","model_eq32_term1_nz",  &model->eq32_term1_nz,    1);

  val = model->A_factored ? 1 : 0;
  H5LTset_attribute_int(   fid,"attrs","model_A_factored",     &val,                     1);
  val = model->try_factor ? 1 : 0;
  H5LTset_attribute_int(   fid,"attrs","model_try_factor",     &val,                     1);
  
// Don't really need all of these (like K), because they're spec'd at CLI

  debug("CHECKP 3: %i\n",(int) status);
  status = attach_double_matrix("gammas", root, gammas);
  status = attach_double_matrix("lhoods", root, lhoods);
  
  if (holdout) {
    debug("CHECKP 3b: %i\n",(int) status);
    H5LTset_attribute_double(fid,"attrs","heldout_bound", &heldout_bound, 1);
    status = attach_double_matrix("heldout_gammas", root, heldout_gammas);
    status = attach_double_matrix("heldout_lhoods", root, heldout_lhoods);
  }

  // The equilibrated factorization of eq32_term1
  // Only built after the first Big Solve and used each subsequent solve
  if (model->A_factored) {
    debug("CHECKP 4a: %i\n",(int) status);
    H5LTset_attribute_int(fid,"attrs","model_A_factor_n", &model->A_factor->n, 1);
    gsl_vector *r    = gsl_vector_calloc(model->A_factor->n);
    gsl_vector *c    = gsl_vector_calloc(model->A_factor->n);
    gsl_vector *ipiv = gsl_vector_calloc(model->A_factor->n);

    for (int i=0; i < model->A_factor->n; i++) { 
      gsl_vector_set(r,   i, model->A_factor->r[i]);
      gsl_vector_set(c,   i, model->A_factor->c[i]);
      gsl_vector_set(ipiv,i, model->A_factor->ipiv[i]);
    }

    status = attach_double_matrix("model_A_factor_af",   root, model->A_factor->af);
    status = attach_double_vector("model_A_factor_ipiv", root, ipiv);
    status = attach_double_vector("model_A_factor_r",    root, r);
    status = attach_double_vector("model_A_factor_c",    root, c);
  }
  
  debug("CHECKP 4b: %i\n",(int) status);

  /* 
     Currently writing eq32_term1 happens after it's built and
     if we don't have HDF5, then we're not checkpointing anyway...so we can assume
     it was written /somehow/ elsewhere -- possibly flat binary.
     and because it doesn't change, no need to checkpoint here.
  */
  //  status = attach_double_matrix("model_eq32_term1", root, model->eq32_term1);

  status = attach_double_matrix("model_mu",         root, model->mu);
  status = attach_double_vector("model_alpha",      root, model->alpha);
  
  debug("CHECKP 5: %i\n",(int) status);
  for (int t=0; t<model->nseq; t++) {
    debug("CHECKP T%i: %i\n", (int) t, (int) status);
    sprintf(roott, "%s/t%i-", root, (int) t);
    
    status = attach_double_matrix("model_influence_doc_weights", roott, model->influence->doc_weights[t]);
    status = attach_double_matrix("model_influence_renormalized_doc_weights",
				  roott, model->influence->renormalized_doc_weights[t]);
  }

  for (int k=0; k<K; k++) {
    debug("CHECKP K%i: %i\n", (int) k, (int) status);
    sprintf(rootk, "%s/k%i-", root, (int) k);

    int    attr_W = 0, attr_T = 0;
    double attr_obs_variance = 0, attr_chain_variance = 0;

    if (model->topic[k]->W != NULL)
      attr_W = model->topic[k]->W;
    if (model->topic[k]->T != NULL)
      attr_T = model->topic[k]->T;
    if (model->topic[k]->obs_variance != NULL)
      attr_obs_variance = model->topic[k]->obs_variance;
    if (model->topic[k]->chain_variance != NULL)
      attr_chain_variance = model->topic[k]->chain_variance;

    sprintf(descriptor, "model_topic_W-k%i", (int) k);
    H5LTset_attribute_int(fid, "attrs",    descriptor, &attr_W, 1);
    sprintf(descriptor, "model_topic_T-k%i", (int) k);
    H5LTset_attribute_int(fid, "attrs",    descriptor, &attr_T, 1);
    sprintf(descriptor, "model_topic_obs_variance-k%i", (int) k);
    H5LTset_attribute_double(fid, "attrs", descriptor, &attr_obs_variance, 1);
    sprintf(descriptor, "model_topic_chain_variance-k%i", (int) k);
    H5LTset_attribute_double(fid, "attrs", descriptor, &attr_chain_variance, 1);

    /**********
    /* These goobers are free()'d after every M step.
    /* They arent NULLed by gsl_*_free(), but they don't need to be stored
    //    status = attach_double_vector("model_topic_zeta", rootk, model->topic[k]->zeta);
    //    status = attach_double_matrix("model_topic_fwd_variance", rootk, model->topic[k]->fwd_variance);
    //    status = attach_double_matrix("model_topic_fwd_mean", rootk, model->topic[k]->fwd_mean);
    //    status = attach_double_vector("model_topic_T_vct", rootk, model->topic[k]->T_vct);
    **********/
    status = attach_double_matrix("model_topic_obs",               rootk, model->topic[k]->obs);
    status = attach_double_matrix("model_topic_e_log_prob",        rootk, model->topic[k]->e_log_prob);
    status = attach_double_matrix("model_topic_mean",              rootk, model->topic[k]->mean);
    status = attach_double_matrix("model_topic_variance",          rootk, model->topic[k]->variance);
    //    status = attach_double_matrix("model_topic_mean_t",            rootk, model->topic[k]->mean_t);
    //    status = attach_double_matrix("model_topic_variance_t",        rootk, model->topic[k]->variance_t);
    status = attach_double_matrix("model_topic_w_phi_l",           rootk, model->topic[k]->w_phi_l);
    status = attach_double_matrix("model_topic_w_phi_sum",         rootk, model->topic[k]->w_phi_sum);
    status = attach_double_matrix("model_topic_w_phi_l_sq",        rootk, model->topic[k]->w_phi_l_sq);
    status = attach_double_matrix("model_topic_m_update_coeff",    rootk, model->topic[k]->m_update_coeff);
    status = attach_double_matrix("model_topic_m_update_coeff_g",  rootk, model->topic[k]->m_update_coeff_g);
    status = attach_double_matrix("model_topic_influence_sum_lgl", rootk, model->influence_sum_lgl[k]);
    status = attach_double_matrix("topic_suffstats",               rootk, topic_suffstats[k]);
  }
 
  debug("CHECKP 7: %i\n", (int) status);
  status = H5Fclose(fid);

  return (int) status;
}
