/*********************************************************************
 Copyright 2007 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

 $Id: nipals.h,v 1.4 2007/05/18 21:55:27 neuhausi Exp $
**********************************************************************/

#ifndef nipals_h
#define nipals_h 

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
EXTERN int verbose;
EXTERN int debug;
#include "files.h"
#include "utils.h"
#include "stats.h"
#include "matrix.h"
#include "methods.h"

void 
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **factor_file,
		char **output_file,
		int *samples,
		int *variables,
		int *components,
		int *tolerance,
		int *all,
		int *transpose);

RESULT **
run_nipals(MATRIX *matrix,
	   int components,
	   int tolerance,
	   int all);

double
is_continue (double *T_old,
	     double **T,
	     int m,
	     int idx);

RESULT **
set_res_object(MATRIX *matrix,
	       int c,
	       int e);

char **
set_out_files(char *output_file,
	      int nres);

#endif

