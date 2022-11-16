/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: files.h,v 1.16 2008/06/11 18:41:38 neuhausi Exp $
**********************************************************************/

#ifndef files_h
#define files_h

#include <malloc.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
EXTERN int verbose;
EXTERN int debug;
#include "utils.h"
#include "stats.h"

typedef struct factor {
  int n;          // number of factors
  int s;          // total number of samples
  int *classes;   // number of levels in each factor 
  int **members;  // class id (i.e, 0, 1, etc) of member in each level in each factor
  int **factors;  // number of members in each factor in each level
  char **samples; // sample names
  char **fnames;  // factor names
  char ***cnames; // levels names
} FACTOR;

typedef struct covariate {
  int n;            // number of covariates
  int s;            // number of samples
  double **values;  // value for the covariate
  char **samples;   // sample names
  char **cnames;    // covariate names
} COVARIATE;

typedef struct group {
  int n;             // number of groups
  int c;             // covariate number id
  int f;             // factor number id
  int *levels;       // number of levels of each group
  int *samples;      // number of samples of each group
  int **groups;      // indeces for the samples in the covariate or factor object
  double **values;  // y values from the covariate object
  char **names;      // level names
} GROUP;

typedef struct vector {
  int n;
  double *data;
} VECTOR;

typedef struct var_vector {
  int n;
  double *mean;
  double *stdev;
} VAR_VECTOR;

typedef struct xdata {
  int n;
  int s;
  int *md;
  char **variables;
  char **samples;
  VECTOR **vectors;
} XDATA;

typedef struct ksdata {
  int w;
  int f;
  int *c;
  char **wells;
  char **features;
  double ***data;  
} KSDATA;

typedef struct matrix {
  int n;
  int m;
  int *md;
  char **rows;
  char **cols;
  double **data;
} MATRIX;

typedef struct nmatrix {
  int n;
  MATRIX **m;
} NMATRIX;

typedef struct result {
  int n;
  int a;
  char **attributes;
  char **variables;
  VECTOR **vectors;
} RESULT;

typedef struct pca_result {
  RESULT *gres;
  RESULT *sres;
} PCA_RESULT;

int
is_binary_file(char *file_name);

int
file_exists(char *filename);

XDATA *
read_xdata(char *file_name, 
	   int n_samples,
	   int n_variables);

XDATA *
read_binary_xdata(char *file_name);

MATRIX *
read_matrix(char *file_name, 
	    int n,
	    int m);

MATRIX *
read_binary_matrix(char *file_name);

FACTOR *
read_factor(char *file_name, 
	    int n_factors,
            int n_samples);

COVARIATE *
read_covariate(char *file_name, 
	       int n_covariates,
	       int n_samples);

GROUP *
set_group(FACTOR *factor,
	  COVARIATE *covariate,
	  char *fit_factor,
	  char *agg_factor,
	  char *zero,
	  int log);

VECTOR *
read_pattern(char *file_name, 
	     int n);

char **
read_file_column(char *file_name,
		 int rows,
                 int header);

void
print_results(char *out_file,
	      RESULT *res,
	      int order,
	      int idx);

void
print_mult_results(int n_res,
		   char **out_file,
		   RESULT **res,
		   int order,
		   int idx);

void
print_PCA_results(char *out_file,
		  PCA_RESULT *res,
		  int axis);

int 
count_lines(char *file_name);

int 
count_wells_ks(char *file_name);

int 
count_fields(char *file_name);

int 
count_not_null_fields(char *file_name);

KSDATA *
read_ksdata(char *file_name,
	    int n_features,
            int n_wells);

#endif
