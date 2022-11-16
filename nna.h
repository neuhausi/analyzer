/*********************************************************************
 Copyright 2006 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus and Robert E. Bruccoleri
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Many thanks to Mark Friedrichs, Carlod Rios and John Hinsdale
   for their help during debugging this code.

 $Id: nna.h,v 1.6 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#ifndef nna_h
#define nna_h 

#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <malloc.h>
EXTERN int verbose;
EXTERN int debug;
#include "files.h"
#include "utils.h"
#include "metrics.h"
#include "var_metrics.h"
#include "methods.h"
#include "stats.h"

enum METRIC_ENUM { EUCLIDEAN, MINKOWSKI, MANHATTAN, MAXIMUM, PEARSON, SPEARMAN, LSPEARMAN };

void
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **profile_file,
		char **output_file,
		int *samples,
		int *variables,
		enum METRIC_ENUM *metrics,
		double *order,
		char ***obj_names,
		int *transpose,
		int *inverse,
		int *n_objects);

VAR_VECTOR *
create_var_pattern(int n_objects, int *obj_numbers, XDATA *xdata);

int
search_obj(XDATA *xdata,
	   char *obj_name);

void
calculate_nna(XDATA *xdata,
	      VECTOR *pattern,
	      char *name,
	      char *output_file,
	      enum METRIC_ENUM metrics,
	      double power,
              int inverse);

void
calculate_nna_var(XDATA *xdata,
		  VAR_VECTOR *var_pattern,
		  char *name,
		  char *output_file,
		  enum METRIC_ENUM metrics,
		  double power,
		  int inverse);

void
print_nna_results(XDATA *xdata,
		  char *name,
		  char *output_file,
		  enum METRIC_ENUM metrics,
		  double sum,
		  double sum2,
		  double *nna);

#endif

