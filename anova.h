/*********************************************************************
 Copyright 2006 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: anova.h,v 1.4 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#ifndef anova_h
#define anova_h 

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
EXTERN int verbose;
EXTERN int debug;
#include "files.h"
#include "utils.h"
#include "stats.h"

enum TT { BETWEEN, WITHIN, MIXED };

void
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **factor_file,
		char **output_file,
		int *samples,
		int *variables,
		int *inter,
		enum TT *tt,
		int *permuts);

RESULT *
run_one_way_anova(XDATA *xdata,
		  FACTOR *factor, 
		  int permuts);

RESULT *
run_multi_factorial_anova(XDATA *xdata,
			  FACTOR *factor, 
			  enum TT tt, 
			  int inter,
			  int permuts);

void
anova_one_way(double data[], int factor[], int n_tot, int cls,
	      double *prob, char **name, double *cls_avg);

void
anova_one_way_within(double data[], int **factor, int n_tot, int cls[], double *prob,
		     char **nameS, char **name);

void
anova_two_way(double data[], int **factor, int n_tot, int cls[],
	      int inter, double *probA, double *probB, double *probAB,
	      char **nameA, char **nameB);

void
anova_two_way_mixed(double data[], int **factor, int n, int cls[],
		    double *probA, double *probB, double *probAB,
		    char **nameS, char **nameA, char **nameB);

void
anova_two_way_within(double data[], int **factor, int n, int cls[],
		     double *probA, double *probB, double *probAB,
		     char **nameS, char **nameA, char **nameB);

void
anova_three_way(double data[], int **factor, int n_tot, int cls[], int inter, 
		double *probA, double *probB, double *probC, 
		double *probAB, double *probAC, double *probBC, double *probABC,
		char **nameA, char **nameB, char **nameC);

#endif

