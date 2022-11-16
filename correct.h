/*********************************************************************
 Copyright 2006 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: correct.h,v 1.3 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#ifndef correct_h
#define correct_h 

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
EXTERN int verbose;
EXTERN int debug;
#include "files.h"
#include "utils.h"
#include "stats.h"
#include "methods.h"

enum TT { FDR, BONFERRONI };

void
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **output_file,
		int *cols,
		int *variables,
		int **inds,
		double *Q,
		enum TT *tt,
		int *n_inds);

RESULT *
correct(MATRIX *matrix,
	double Q,
	enum TT tt, 
	int *inds,
	int n_cols,
	double **cutoff,
	int  **nsig);

void
print_FDR_summary(char *out_file,
		  RESULT *res,
		  double *cutoff,
		  int *nsig,
		  int n_cols);

#endif

