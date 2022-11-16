/*********************************************************************
 Copyright 2006 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: svd.h,v 1.3 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#ifndef svd_h
#define svd_h 

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

enum METHOD_ENUM { COR, COV };

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
		int *transpose);

RESULT *
run_svd(MATRIX *matrix,
	int components,
	int eigenvs,
	enum METHOD_ENUM method);

RESULT *
set_res_object(int s,
	       int n,
	       int c,
	       int e);

void
corcol (double **data, int n, int m);

void
covcol (double **data, int n, int m);

#endif

