/*********************************************************************
 Copyright 2004 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

 $Id: mixed.h,v 1.17 2010/04/09 22:07:28 neuhausi Exp $
**********************************************************************/

#ifndef mixed_h
#define mixed_h 

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

enum METHOD { REML, ML, MIVQUE };

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
  int **rzeros;      // valid columns in random matrix
  int **rvalid;      // valid lsm values in random matrix
  int **fzeros;      // valid columns in fixed matrix
  int **fvalid;      // valid lsm values in fixed matrix

  char *formula;     // formula

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
		enum METHOD *method,
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

MATRIX *
set_G_derivatives(MATRIX *G,
		  int **indeces,
		  int index);

void
update_G_matrix(MATRIX **G,
		MATRIX *parameters,
		int **indeces);

void
update_R_matrix(MATRIX **R,
		MATRIX *parameters);

RESULT *
calculate_model(XDATA *xdata,
		FACTOR *factor,
		COVARIATE *covariate,
		char **variable_names,
		int variable_number,
		MODEL *model,
		enum METHOD method,
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


MATRIX *
set_coef_matrix(MATRIX *XRX,
		MATRIX *XRZ,
		MATRIX *ZRX,
		MATRIX *ZRZG);

MATRIX *
set_initial_W_matrix(MATRIX *X,
		     MATRIX *Z,
		     MATRIX *Y,
		     MATRIX **bu,
		     MATRIX **varb);

MATRIX *
compute_W_transform(MATRIX *W0a,
		    MATRIX *G,
		    MATRIX *X,
		    MATRIX *R,
		    MATRIX **bu,
		    MATRIX **varb,
		    enum METHOD method,
		    double *ml,
		    double *reml,
		    int n);

MATRIX *
compute_W0_augment(MATRIX *W0,
		   MATRIX *iR,
		   MATRIX *X,
		   MATRIX *Z,
		   MATRIX *G);

MATRIX *
compute_W0_transform(MATRIX *X,
		     MATRIX *Z,
		     MATRIX *iR,
		     MATRIX *Y);

MATRIX *
set_W_trans(MATRIX *X,
	    MATRIX *iR,
	    MATRIX *Y);

MATRIX *
set_y_border(MATRIX *XRY,
	     MATRIX *ZRY);


MATRIX *
compute_subvariance_matrix(MATRIX *iR,
			   MATRIX *Z,
			   MATRIX *G);

MATRIX *
compute_variance_matrix(MATRIX *iR,
			MATRIX *Z,
			MATRIX *M);

MATRIX *
compute_X_star_matrix(MATRIX *V);
		     
void
compute_derivatives(MATRIX *X,
		    MATRIX *Z,
		    MATRIX *iR,
		    MATRIX *G,
		    MATRIX *Y,
		    MATRIX **hess,
		    MATRIX **grad,
		    MATRIX *W,
		    MATRIX *bu,
		    MATRIX *varb,
		    NMATRIX *G_deriv,
		    NMATRIX *R_deriv,
		    enum METHOD method,
		    int niter);

MATRIX *
get_Newton_direction_matrix(MATRIX **hess,
			    MATRIX *grad,
			    int p);

void
Newton_Raphson(MATRIX *X,
	       MATRIX *Z,
	       MATRIX *iR,
	       MATRIX *G,
	       MATRIX *R,
	       int **rindeces,
	       MATRIX *Y,
	       MATRIX *parameters,
	       MATRIX **hess,
	       MATRIX **grad,
	       MATRIX **bu,
	       MATRIX **varb,
	       NMATRIX *G_deriv,
	       NMATRIX *R_deriv,
	       enum METHOD method,
	       double *ml,
	       double *reml);

MATRIX * 
compute_EM_matrix(MATRIX *X,
	          MATRIX *Z,
	          MATRIX *G,
	          MATRIX *Y,
		  MATRIX *parameters);
		  	      
void
Expected_Maxim(MATRIX *X,
	       MATRIX *Z,
	       MATRIX *iR,
	       MATRIX *G,
	       MATRIX *R,
	       int **rindeces,
	       MATRIX *Y,
	       MATRIX *parameters,
	       MATRIX **hess,
	       MATRIX **grad,
	       MATRIX **bu,
	       MATRIX **varb,
	       NMATRIX *G_deriv,
	       NMATRIX *R_deriv,
	       enum METHOD method,
	       double *ml,
	       double *reml);

NMATRIX *
set_contrast_matrices(MATRIX *des,
		      NMATRIX *hyp,
                      int type);

MATRIX *
sas_ortogonalize(MATRIX *ef,
		 MATRIX *in);

NMATRIX *
set_ef3_matrix(MODEL *model,
	       MATRIX *gef,
               int *rank,
               int **indeces,
	       int skip);

MATRIX *
set_gef_matrix(MATRIX *x,
	       int **rank);

MATRIX *
set_design_matrix(FACTOR *factor, 
		  COVARIATE *covariate,
		  MODEL *model,
                  int ***indeces);

MATRIX *
set_random_matrix(FACTOR *factor,
		  COVARIATE *covariate,
		  MODEL *model,
		  int ***rindeces);

MATRIX *
set_fixed_matrix(FACTOR *factor,
		 COVARIATE *covariate,
		 MODEL *model,
		 int ***rindeces);

void
mixed(VECTOR *y,
      MATRIX *z,
      int **rindeces,
      MATRIX *f,
      MATRIX *R,
      MATRIX *G,
      MATRIX **hess,
      MATRIX *grad,
      NMATRIX *G_deriv,
      NMATRIX *R_deriv,
      MATRIX *design,
      int **indeces,
      MATRIX *X,
      NMATRIX *H,
      int *ndf,
      NMATRIX *C,
      NMATRIX *lsm,
      RESULT *res,
      MODEL *model,
      enum METHOD method,
      int idx,
      char *reference,
      int lsmo,
      int islog);

float
compute_f_statistic(MATRIX *H,
		    MATRIX *varb,
		    MATRIX *b);

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

int *
calculate_dof(MATRIX *des,
	      MODEL *model,
	      int **indeces,
	      int ssq);

int *
get_contained_dof(MATRIX *des,
		  MODEL *model,
		  int **indeces,
		  int *ndf);

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
