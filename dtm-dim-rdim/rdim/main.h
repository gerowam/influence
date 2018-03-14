#ifndef MAINH
#define MAINH

#include <stdlib.h>
#include <string.h>
#include "data.h"
#include "lda-seq.h"
#include "lda.h"
#include <gsl/gsl_matrix.h>

#ifdef CUDA
#include <cuda_runtime.h>
#include <cublas.h>
#include <cusparse.h>

extern cusparseHandle_t hndl;
#endif

typedef struct dtm_fit_params {
    char* datafile;
    char* outname;
    char* heldout;
    int start;
    int end;
    int ntopics;
    int lda_max_em_iter;
    double top_obs_var;
    double top_chain_var;
    double alpha;
} dtm_fit_params;

void set_thread_affinity();

#endif
