/*********************************************************************
 Copyright 2006 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: ttest.h,v 1.3 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#ifndef ttest_h
#define ttest_h 

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
EXTERN int verbose;
EXTERN int debug;
#include "files.h"
#include "utils.h"
#include "stats.h"

enum TT { PAIRED, UNEQUAL, EQUAL };

void
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **factor_file,
		char **output_file,
		int *samples,
		int *variables,
		enum TT *tt,
		double *md,
		int *permuts);

RESULT *
run_ttest(XDATA *xdata,
	  FACTOR *factor,
	  enum TT tt, 
	  double md,
	  int permuts);

void
ttest(double data1[], int n1, double data2[], int n2,
      double *t, double *df, double *prob,
      double *m1, double *v1, double *m2, double *v2);

void
tutest(double data1[], int n1, double data2[], int n2,
       double *t, double *df, double *prob,
       double *m1, double *v1, double *m2, double *v2);

void
tptest(double data1[], double data2[], int n,
       double *t, double *df, double *prob,
       double *m1, double *v1, double *m2, double *v2);

#endif

