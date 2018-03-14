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

#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "hdf5.h"
#include "hdf5_hl.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "param.h"
#include "gsl-wrappers.h"
#include "lda-seq.h"

gsl_matrix* read_double_matrix(const char* descriptor, const char* dir);
gsl_vector* read_double_vector(const char* descriptor, const char* dir);

int read_checkpoint(const char*   indir,            lda_seq *&model,
		    gsl_matrix**  topic_suffstats,  gsl_matrix *&gammas,
		    gsl_matrix   *&lhoods,           gsl_matrix *&heldout_gammas,
		    gsl_matrix   *&heldout_lhoods,   double&     bound,
		    double&       heldout_bound,    double&     convergence,
		    double&       old_bound,        int&        iter,
		    short&        final_iters_flag, const bool  holdout,
 	            unsigned int& last_iter,        const int   K);

herr_t attach_double_matrix(const char* descriptor, const char* dir, gsl_matrix* mat);
herr_t attach_double_vector(const char* descriptor, const char* dir, gsl_vector* vec);

int write_checkpoint(const char*        outdir,           lda_seq*     model,
		     gsl_matrix**       topic_suffstats,  gsl_matrix*  gammas,
		     gsl_matrix*        lhoods,           gsl_matrix*  heldout_gammas,
		     gsl_matrix*        heldout_lhoods,   const double bound,
		     const double       heldout_bound,    const double convergence,
		     const double       old_bound,        const int    iter,
		     const short        final_iters_flag, const bool   holdout,
		     const unsigned int last_iter,        const int    K);

#endif
