/*********************************************************************
 Copyright 2006 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

 $Id: fit.h,v 1.13 2007/07/19 22:17:23 neuhausi Exp $
**********************************************************************/

#ifndef fit_h
#define fit_h 

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
EXTERN int verbose;
EXTERN int debug;
#include "files.h"
#include "utils.h"
#include "methods.h"
#include "stats.h"
#include "solution.h"
#include "function.h"

enum FUNCS { NLLS4P };

void 
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **factor_file,
		char **covariate_file,
		char **output_file,
		int *samples,
		int *variables,
		char ***variable_names,
		int *variable_number,
		char **fit_factor,
		char **agg_factor,
		char **zero,
		int *logtrans,
		enum FUNCS *func,
		double **params,
		double *md,
		int *progress);

RESULT *
set_result(GROUP *group,
	   int n);

void
set_attribute_result(int i,
		     char *str,
		     RESULT *res);

RESULT *
fit_data(XDATA *xdata,
	 FACTOR *factor,
	 COVARIATE *covariate,
	 char **variable_names,
	 int variable_number,
	 char *fit_factor,
	 char *agg_factor,
	 char *zero,
	 int logtrans,
	 enum FUNCS func,
	 double *params,
	 double md,
	 int progress);

void
no_fit(VECTOR *vr,
       int gr,
       char *variable_name);

int
search_variable(XDATA *xdata,
		char *variable_name);

void
fit_nlsr4p(VECTOR *v,
	   GROUP *group,
	   double *params,
	   RESULT *res,
	   int idx,
	   int gr, 
	   int progress);

#endif
