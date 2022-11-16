/*********************************************************************
 Copyright 2006 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Principal Components Analysis or the Karhunen-Loeve expansion is a
   classical method for dimensionality reduction or exploratory data
   analysis.  This code was adapted from a public domain version by
   F. Murtagh and A. Heck, Multivariate Data Analysis, Kluwer
   Academic, Dordrecht, 1987 . Also several subrutines were also
   adapted from Numerical Recipes in C, W H Press et al., 1988

 $Id: pca.h,v 1.4 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#ifndef pca_h
#define pca_h 

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
EXTERN int verbose;
EXTERN int debug;
#include "files.h"
#include "utils.h"
#include "stats.h"
#include "solution.h"
#include "methods.h"

enum METHOD_ENUM { COR, COV, SSCP };
enum AXIS { VARIABLES, SAMPLES, BOTH };

void 
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **factor_file,
		char **output_file,
		int *samples,
		int *variables,
		int *components,
		int *eigenvs,
		enum METHOD_ENUM *method,
		enum AXIS *axis,
		int *transpose);

PCA_RESULT *
run_pca(MATRIX *matrix,
	int components,
	int eigenvs,
	enum METHOD_ENUM method,
	enum AXIS axis);

PCA_RESULT *
set_res_object(int s,
	       int n,
	       int c,
	       int e);

void
corcol (double **data, int n, int m, double **symmat);

void
covcol (double **data, int n, int m, double **symmat);

void
scpcol (double **data, int n, int m, double **symmat);

#endif

