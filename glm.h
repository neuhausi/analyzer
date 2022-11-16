/*********************************************************************
 Copyright 2004 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

 $Id: glm.h,v 1.16 2009/02/11 22:40:20 neuhausi Exp $
**********************************************************************/

#ifndef glm_h
#define glm_h 

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
EXTERN int verbose;
EXTERN int debug;
#include "files.h"
#include "utils.h"
#include "methods.h"
#include "stats.h"
#include "matrix.h"
#include "solution.h"

typedef struct model {

  int *infactor;     // factors included in the model
  int *infactorrand; // random factors included in the model 
  int *infactornest; // nested factors included in the model 

  int n;             // number of nested factors
  int **nested;      // nested factors id(s)
  int *numnest;      // number of factors in the nested factor
  int *denomnest;    // denominator for the nested factor

  int *incovariate;  // covariates included in the model

  int f;             // number of fixed effects
  int *fixed;        // fixed effect id (fac or cov)
  int *tyfixed;      // type of fixed effect (fac=0, cov=1)

  int i;             // number of interactions
  int **inter;       // term id (fac or cov)
  int **tyinter;     // type of term (fac=0, cov=1)
  int *numinter;     // number of terms in each interaction
  int *maxnumlev;    // maximum number of levels in each interaction

  int r;             // number of random effects
  int *rand;         // factor number
  int *randinter;    // random interations

  int s;             // number samples
  int k;             // number of samples skipped due to Unassigned factors or covariates
  int *include;      // Flag to indicate inclusion or exclussion of samples (0 or 1)

  int t;             // number of terms (fixed + inter)
  int *numdep;       // number of dependencies
  int **dep;         // term id for dependency
  int *randdep;      // random term

  char **factors;    // factors in the model
  int ***lsm;        // ids for the lsm matrices
  char ***levels;    // levels in the model
  int *termlen;      // number of columns in each term
  int **zeros;       // valid columns
  int **valid;       // valid lsm values

  char *formula;     // fornula

} MODEL;

void 
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **factor_file,
		char **covariate_file,
		char **output_file,
		int *samples,
		int *variables,
		char ***variable_name,
		int *variable_number,
		char ***fixed_effs,
		int *n_fixed_effs,
		char ***interacts,
		int *n_interacts,
		char ***random_facs,
		int *n_random_facs,
		int *ssq,
		double *Q,
		char **reference,
		int *lsmo,
		int *islog,
		int *permuts);

RESULT *
set_result(int ndata,
	   FACTOR *factor,
	   COVARIATE *covariate,
	   MODEL *model,
	   MATRIX *des,
	   double Q,
	   char *reference,
	   int lsmo);

void
set_attribute_result(int i,
		     char *str,
		     RESULT *res);

void
validate_reference(MODEL *model,
		   char *reference);

void
mult_test_correction(RESULT *res,
		     MODEL *model,
		     double Q);

RESULT *
calculate_glm(XDATA *xdata,
	      FACTOR *factor,
	      COVARIATE *covariate,
	      char **variable_names,
	      int variable_number,
	      MODEL *model,
	      int ssq,
	      double Q,
	      char *reference,
	      int lsmo,
	      int islog,
	      int permuts);

int
search_variable(XDATA *xdata,
		char *variable_name);

MATRIX *
set_X_matrix(MATRIX *des);

NMATRIX *
set_contrast_matrices(MATRIX *des,
		      NMATRIX *hyp);

MATRIX *
sas_ortogonalize(MATRIX *ef,
		 MATRIX *in);

NMATRIX *
set_ef3_matrix(MODEL *model,
	       MATRIX *gef,
               int *rank,
               int **indeces);

MATRIX *
set_gef_matrix(MATRIX *x,
	       int **rank);

MATRIX *
set_design_matrix(FACTOR *factor, 
		  COVARIATE *covariate,
		  MODEL *model,
                  int ***indeces);

MATRIX *
clean_design_matrix(int *rna,
		    MATRIX *des,
		    int na);

NMATRIX *
set_lsmeans_matrix(MODEL *model,
		   MATRIX *des,
		   int **indeces);

MATRIX *
design(FACTOR *factor,
       COVARIATE *covariate,
       MODEL *model,
       int o,
       int t);

int
is_nested(FACTOR *factor,
	  COVARIATE *covariate,
	  MODEL *m,
	  int o);

MODEL *
set_model(FACTOR *factor,
	  COVARIATE *covariate,
	  char **fixed_effs,
	  int n_fixed_effs,
	  char **interacts,
	  int n_interacts,
	  char **random_facs,
	  int n_random_facs);

double
get_sum_of_squares(MATRIX *H,
		   MATRIX *C,
		   MATRIX *b);

MATRIX *
get_least_square_means(MATRIX *L,
		       MATRIX *B);

void
parameter_estimates(VECTOR *Y,
		    MATRIX *des,
		    MATRIX *X,
		    MATRIX **b,
		    double *sse);

void
satterthwaite_approx (MODEL *model,
		      double **temp,
		      double *df,
		      double *f,
		      double *ms,
		      int idx);

void
glm(VECTOR *y,
    MATRIX *design,
    MATRIX *X,
    NMATRIX *H,
    NMATRIX *C,
    NMATRIX *lsm,
    RESULT *res,
    MODEL *model,
    int idx,
    char *reference,
    int lsmo,
    int islog);

#endif
