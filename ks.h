/*********************************************************************
 Copyright 2006 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: ks.h,v 1.4 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#ifndef ks_h
#define ks_h 

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
EXTERN int verbose;
EXTERN int debug;
#include "files.h"
#include "utils.h"
#include "stats.h"

void
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **output_file,
		int *wells,
		int *features,
		char ***reference,
		int *n_reference,
		char **feature_name);

void
kstwo(double data1[], int n1, double data2[], int n2, double *d);

RESULT *
set_result(KSDATA *ksdata,
	   char *feature_name);

RESULT *
calculate_ks (KSDATA *ksdata,
	      char **reference,
	      int n_reference,
	      char *feature_name);

#endif
