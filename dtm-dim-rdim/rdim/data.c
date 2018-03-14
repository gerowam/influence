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

#include <string.h>
#include "gflags.h"
#include "data.h"

// Rounding precision when reading in ASCII data
#define PREC 15

DEFINE_double(sigma_l,    0.05, "DIM / rDIM hyper-param (see Gerrish & Blei).");
DEFINE_double(sigma_d,    0.05, "DTM hyper-param (see Blei & Lafferty).");
DEFINE_double(sigma_c,    0.05, "DTM hyper-param (see Blei & Lafferty).");
DEFINE_double(sigma_cv,   1e-6, "DTM hyper-param (see Blei & Lafferty).");
DEFINE_double(resolution, 1,    "The resolution. Used to determine how far out the beta mean should be.");
DEFINE_int32(max_number_time_points, 200,  "Used for the influence window.");
DEFINE_double(time_resolution,       0.5,  "This is the number of years per time slice.");
DEFINE_double(influence_mean_years,  20.0, "How many years is the mean number of citations?");
DEFINE_double(influence_stdev_years, 15.0, "How many years is the stdev number of citations?");
DEFINE_int32(influence_flat_years,  -1,    "How many years is the influence nonzero? If <=0, a lognormal is used.");
DECLARE_string(normalize_docs);
DECLARE_string(model);
DECLARE_int32(debug);

/*
 * seq corpus range: [start, end)
 * creates a subset of time slices
 *
 */
corpus_seq_t* make_corpus_seq_subset(corpus_seq_t* all, int start, int end) {
    int n;
    corpus_seq_t* subset_corpus = (corpus_seq_t*) malloc(sizeof(corpus_seq_t));

    subset_corpus->nterms = all->nterms;
    subset_corpus->len    = end - start;
    subset_corpus->ndocs  = 0;
    subset_corpus->corpus = (corpus_t**) malloc(sizeof(corpus_t*) * subset_corpus->len);

    if (all->nmetafields <= 0) {
        outlog("%s","FATAL: Because this is rDIM, you need to specify --metafields as a positive integer\n");
        exit(1);
    }

    subset_corpus->nmetafields = all->nmetafields; // Added or rDIM by AG

    for (n = start; n < end; n++) {
        subset_corpus->corpus[n - start] = all->corpus[n];
        subset_corpus->ndocs            += all->corpus[n]->ndocs;
    }

    return(subset_corpus);
}

/*
 * collapse a sequential corpus to a flat corpus
 *
 */
corpus_t* collapse_corpus_seq(corpus_seq_t* c) {
  corpus_t* collapsed = (corpus_t*) malloc(sizeof(corpus_t));
  collapsed->ndocs  = c->ndocs;
  collapsed->nterms = c->nterms;
  collapsed->doc    = (doc_t**) malloc(sizeof(doc_t*) * c->ndocs);
  collapsed->max_unique = 0;
  int t, n, doc_idx = 0;

  for (t = 0; t < c->len; t++) {
    for (n = 0; n < c->corpus[t]->ndocs; n++) {
      collapsed->doc[doc_idx] = c->corpus[t]->doc[n];

      if (collapsed->doc[doc_idx]->nterms > collapsed->max_unique)
	collapsed->max_unique = collapsed->doc[doc_idx]->nterms;

      doc_idx++;
    }
  }

  assert(doc_idx == collapsed->ndocs);
  return(collapsed);
}

/*
 * read corpus
 *
 */
corpus_t* read_corpus(const char* name) {
  int length, count, word, n;
  corpus_t* c;
  char filename[400];
  sprintf(filename, "%s-mult.dat", name);
  outlog("Reading corpus from %s\n", filename);

  c = (corpus_t*) malloc(sizeof(corpus_t));
  c->max_unique = 0;
  FILE* fileptr = fopen(filename, "r");

  if (fileptr == NULL) {
    outlog("FATAL: reading corpus prefix %s. Failing.\n", filename);
    exit(1);
  }

  c->ndocs  = 0; 
  c->nterms = 0;
  c->doc = (doc_t**) malloc(sizeof(doc_t*));
  int grand_total = 0;

  while ((fscanf(fileptr, "%10d", &length) != EOF)) {
    if (length > c->max_unique) 
      c->max_unique = length;

    c->doc = (doc_t**) realloc(c->doc, sizeof(doc_t*)*(c->ndocs+1));
    c->doc[c->ndocs] = (doc_t*) malloc(sizeof(doc_t));
    c->doc[c->ndocs]->nterms = length;
    c->doc[c->ndocs]->total  = 0;
    c->doc[c->ndocs]->log_likelihood = 0.0;

    c->doc[c->ndocs]->word   = (int*)    malloc(sizeof(int)   *length);
    c->doc[c->ndocs]->count  = (int*)    malloc(sizeof(int)   *length);
    c->doc[c->ndocs]->lambda = (double*) malloc(sizeof(double)*length);
    c->doc[c->ndocs]->log_likelihoods = (double*) malloc(sizeof(double)*length);

    for (n = 0; n < length; n++) {
      fscanf(fileptr, "%10d:%10d", &word, &count);
      word = word - OFFSET;

      if (FLAGS_normalize_docs == "occurrence")
	count = 1;

      c->doc[c->ndocs]->word[n]  = word;
      c->doc[c->ndocs]->count[n] = count;
      c->doc[c->ndocs]->total   += count;

      // Is there a better value for initializing lambda?
      c->doc[c->ndocs]->lambda[n] = 0.0;
      c->doc[c->ndocs]->log_likelihoods[n] = 0.0;

      if (word >= c->nterms) 
	c->nterms = word + 1; 
    }

    grand_total += c->doc[c->ndocs]->total;
    c->ndocs = c->ndocs + 1;
  }

  fclose(fileptr);
  outlog("Succesfully read corpus (ndocs = %d; nterms = %d; nwords = %d)\n", c->ndocs, c->nterms, grand_total);

  return(c);
}

/*
 * read corpus sequence
 *
 */
corpus_seq_t* read_corpus_seq(const char* name, const int ntopics, const int nmetafields) {
  char          filename[400];
  corpus_seq_t* corpus_seq = (corpus_seq_t*) malloc(sizeof(corpus_seq_t));

  // read corpus
  corpus_t* raw_corpus = read_corpus(name);
  corpus_seq->nterms   = raw_corpus->nterms;

  // read sequence information
  sprintf(filename, "%s-seq.dat", name);
  outlog("Reading corpus sequence %s.\n", filename);
  FILE* fileptr = fopen(filename, "r");

  if (!fileptr) {
    outlog("FATAL: Error opening dtm sequence file %s.\n", filename);
    exit(1);
  }

  fscanf(fileptr, "%d", &(corpus_seq->len));
  corpus_seq->corpus = (corpus_t**)   malloc(sizeof(corpus_t*)   * corpus_seq->len);

  // allocate corpora
  int doc_idx = 0;
  int ndocs;
  corpus_seq->ndocs = 0;

  // Added by AG:
  if (FLAGS_model == "rdim") {
    if (nmetafields <= 0) {
      outlog("%s","FATAL: Because this is rDIM, you need to specify --metafields as a positive integer.\n");
      exit(1);
    }

    outlog("Reading rDIM meta-data...\n%s","");

#ifdef SPARSE
    corpus_seq->tau = (gsl_spmatrix**) malloc(sizeof(gsl_spmatrix*) * corpus_seq->len); 
#else
    corpus_seq->tau = (gsl_matrix**) malloc(sizeof(gsl_matrix*) * corpus_seq->len); 
#endif

    corpus_seq->nmetafields = nmetafields;
  }
    
  for (int i = 0; i < corpus_seq->len; ++i) {
    // Somewhere in here we need to set corpus_seq->corpus[i]->nterms
    fscanf(fileptr, "%d", &ndocs);
    corpus_seq->ndocs            += ndocs;
    corpus_seq->corpus[i]         = (corpus_t*) malloc(sizeof(corpus_t));
    corpus_seq->corpus[i]->ndocs  = ndocs;
    corpus_seq->corpus[i]->doc    = (doc_t**) malloc(sizeof(doc_t*) * ndocs);

    // Added by AG
    if (FLAGS_model == "rdim") 
      corpus_seq->tau[i]         = read_tau(name, i, ndocs, corpus_seq->nmetafields); 
	
    for (int j = 0; j < ndocs; j++) {
      if (doc_idx >= raw_corpus->ndocs) {
	outlog("Error: too few documents listed in dtm sequence file %s.\n", filename);
	outlog("Current line: %d %d %d.\n", doc_idx, ndocs, j);
	exit(1);
      }

      corpus_seq->corpus[i]->doc[j] = raw_corpus->doc[doc_idx];
      doc_idx++;
    }

    corpus_seq->corpus[i]->nterms = raw_corpus->nterms;
  }
    
  corpus_seq->max_nterms = compute_max_nterms(corpus_seq);

  if (FLAGS_model == "rdim") {
    outlog("Sucessfully read corpus of (time-slices = %d; ndocs = %d; nterms = %d, metafields = %d)\n", 
	   corpus_seq->len, corpus_seq->ndocs, corpus_seq->nterms, corpus_seq->nmetafields);
  } else {
    outlog("Sucessfully read corpus of (time-slices = %d; ndocs = %d; nterms = %d)\n", 
	   corpus_seq->len, corpus_seq->ndocs, corpus_seq->nterms);
  }
    
  free(raw_corpus->doc); // added by YH
  free(raw_corpus);      // added by YH
    
  return(corpus_seq);
}

/*
 * write sequential corpus
 *
 */
void write_corpus_seq(corpus_seq_t* c, char* name) {
    char tmp_string[400];
    int n;

    outlog("writing %d slices to %s (%d total docs)\n", c->len, name, c->ndocs);
    sprintf(tmp_string, "%s-seq.dat", name);
    FILE* seq_file = fopen(tmp_string, "w");
    fprintf(seq_file, "%d", c->len);

    for (n = 0; n < c->len; n++)
        fprintf(seq_file, " %d", c->corpus[n]->ndocs);

    fclose(seq_file);

    corpus_t* flat = collapse_corpus_seq(c);
    sprintf(tmp_string, "%s-mult.dat", name);
    write_corpus(flat, tmp_string);
}

/*
 * write corpus
 *
 */
void write_corpus(corpus_t* c, char* filename) {
  int i, j;
  FILE * fileptr;
  doc_t * d;
  outlog("writing %d docs to %s\n", c->ndocs, filename);
  fileptr = fopen(filename, "w");

  for (i = 0; i < c->ndocs; i++) {
    d = c->doc[i];
    fprintf(fileptr, "%d", d->nterms);

    for (j = 0; j < d->nterms; j++) 
      fprintf(fileptr, " %d:%d", d->word[j], d->count[j]);

    fprintf(fileptr, "\n");
  }

  fclose(fileptr);
}

/*
 * compute the maximum nterms in a corpus sequence
 *
 */
int compute_max_nterms(const corpus_seq_t* c) {
  int max = 0;

  for (int i = 0; i < c->len; i++) {
    corpus_t* corpus = c->corpus[i];

    for (int j = 0; j < corpus->ndocs; j++)
      if (corpus->doc[j]->nterms > max)
	max = corpus->doc[j]->nterms;
  }

  return(max);
}

/*
 * compute the total matrix of counts (W x T)
 *
 */
gsl_matrix* compute_total_counts(const corpus_seq_t* c) {
  int t, d, n;
  gsl_matrix* ret = gsl_matrix_alloc(c->nterms, c->len);

  for (t = 0; t < c->len; t++) {
    corpus_t* corpus = c->corpus[t];

    for (d = 0; d < corpus->ndocs; d++) {
      doc_t* doc = corpus->doc[d];

      for (n = 0; n < doc->nterms; n++)
	minc(ret, doc->word[n], t, (double) doc->count[n]);
    }
  }

  return(ret);
}

// Creates a new array of doubles with kScaledBetaMax elements.
double* NewScaledInfluence(int size) {
  outlog("%s","Creating new scaled influence vector\n");

  double* scaled_influence = new double[size];
  
  if (FLAGS_influence_flat_years > 0) {
    // Note that we round up, to make sure we have at least one epoch.
    int    number_epochs = FLAGS_influence_flat_years * FLAGS_time_resolution;
    double epoch_weight  = 1.0 / number_epochs;

    for (int i=0; i < number_epochs; ++i) 
      scaled_influence[i] = epoch_weight;

    for (int i=number_epochs; i < size; ++i)
      scaled_influence[i] = 0.0;

    return scaled_influence;
  } 
  else if (FLAGS_influence_mean_years > 0) {  // Handle the log-normal distribution.
    double total = 0.0;

    // Here, we're interested more in the median.
    // So we treat the variable mean as median and note this in our paper.
    double scaled_influence_mean     = FLAGS_influence_mean_years;
    double scaled_influence_variance = (FLAGS_influence_stdev_years * FLAGS_influence_stdev_years);
    double tmp                       = (1.0 + (scaled_influence_variance  / 
					       (scaled_influence_mean * scaled_influence_mean)));
    double lognormal_sigma_squared   = log(tmp);
    double lognormal_mu              = (log(scaled_influence_mean) - 0.5 * lognormal_sigma_squared);

    for (int i = 0; i < size; ++i) {
      // Shift right by half a timeframe to avoid corner cases.
      double x    = (i / FLAGS_time_resolution) + (1.0 / FLAGS_time_resolution) / 2;
      double tmp2 = (log(x) - lognormal_mu);
      scaled_influence[i] = (1.0 / (x * sqrt(lognormal_sigma_squared * 2 * 3.1415926))
			     * exp(-tmp2 * tmp2 / (2.0 * lognormal_sigma_squared)));
      total += scaled_influence[i];
    }
  
    for (int i = 0; i < kScaledInfluenceMax; ++i)
      scaled_influence[i] /= total;

    return scaled_influence;
  } else {
    outlog("%s","FATAL: Influence envolope improperly specified\nUse -influence_flat_years, -influence_mean_years or -influence_stdev_years.\n");
    exit(1);
  }
}

#ifdef SPARSE
double g_round_data(double x, unsigned int digits) {
  double fac = pow(10, digits);
  return round(x*fac) / fac;
}

void print_spmatrix_data(gsl_spmatrix *M) {
  for (int row = 0; row < M->size1; row++)
    for (int col = 0; col < M->size2; col++) {
      outlog("%.14f\n", gsl_spmatrix_get(M,row,col));
    }
}

// Read in a GSL matrix defining the metadata for a given time-slice.
// These are pre-computed by the `compile_metadata.py` script
// Added by AG
gsl_spmatrix* read_tau(const char* prefix, const int t, 
		       const int   ndocs,  const int nmetafields) { // SPARSE
  gsl_spmatrix* ret = gsl_spmatrix_alloc(ndocs, nmetafields); // Number of docs in corpus[t]
  gsl_spmatrix_set_zero(ret);
  char  *line = NULL;
  size_t len  = 0;

  char filename[400];
  sprintf(filename, "/%s-metadata/tau-%i.sparse", prefix, t); 
  FILE* fptr = fopen(filename, "r");
  
  outlog("Reading in metadata (%s), tau, at time %i...", filename, t);
  
  while (getline(&line, &len, fptr) != -1) {
    char *i = strtok(line, " ");
    char *j = strtok(NULL, " ");
    char *v = strtok(NULL, " "); // Not terribly safe, but we're pretty sure these files are well-formed.
    
    // Force the dimensions set by the user (even if there are cells in the DATs that exceed them.)
    if (atoi(j) >= nmetafields) continue;

    gsl_spmatrix_set(ret, (size_t) atoi(i), (size_t) atoi(j), g_round_data(strtod(v, NULL), PREC));
    line = NULL;
  }
  
  fclose(fptr);
  
  gsl_spmatrix *M = gsl_spmatrix_compcol(ret);
  gsl_spmatrix_free(ret);
  
  outlog("Success: (%d x %d)\n", (int) M->size1, (int) M->size2);
  
  return(M);
}

#else
// Read in a GSL matrix defining the metadata for a given time-slice.
// These are pre-computed by the `compile_metadata.py` script
// Added by AG
gsl_matrix* read_tau(const char* prefix, const int t, const int ndocs, const int nmetafields) { // DENSE
  gsl_matrix* M = gsl_matrix_alloc(ndocs, nmetafields); // Number of docs in corpus[t]

  char filename[400];
  sprintf(filename, "/%s-metadata/tau-%i.dat", prefix, t); 
  FILE* fptr = fopen(filename, "r");

  outlog("Reading in metadata (%s), tau, at time %i...", filename, t);    

  if (gsl_matrix_fscanf(fptr, M) == 0) {
    outlog("Success: (%d x %d)\n", (int) M->size1, (int) M->size2);
    fclose(fptr);
  } else {
    outlog("FATAL: failed to read in metadata (%s), tau, at time %i", filename, t);
    fclose(fptr);
    exit(1);
  }

  if (FLAGS_debug) gsl_matrix_fprintf(stderr, M, "%.14f");

  return(M);
}
#endif
