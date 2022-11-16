/*********************************************************************
 Copyright 2006 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: hyper.h,v 1.6 2009/01/09 15:18:46 neuhausi Exp $
**********************************************************************/

#ifndef hyper_h
#define hyper_h 

#include <stdio.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
EXTERN int verbose;
EXTERN int debug;
#include "files.h"
#include "utils.h"
#include "stats.h"

enum TAIL { AUTO, LOWER, UPPER };

void
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **list_file,
		char **output_file,
		char **matrix_file,
		int *cols,
		int *rows,
		int *n_objects,
		char **obj_name,
		int *hits,
		double *Q,
		enum TAIL *tail);

RESULT *
set_result(XDATA *xdata);

int
search_object(XDATA *xdata,
	      char *obj_name);


void
uppercase(char *str,
	  int len);

void
set_vector(XDATA *xdata,
	   char **objects,
	   int n_objects,
	   VECTOR *v);

RESULT *
calculate_hyper (XDATA *xdata,
		 char *obj_name,
		 char **objects,
		 int n_objects,
		 int *hits,
		 double Q,
		 enum TAIL tail);

void
hyper (double q,
       double m,
       double n,
       double k,
       RESULT *res,
       int i,
       enum TAIL tail);

double
dhyper (double q,
	double m,
	double n,
	double k);

double
pdhyper (double q,
	 double m,
	 double n,
	 double k);


double
dbinom_raw (double q, 
           double k, 
           double p, 
           double pq);

double 
stirlerr (double n);

double 
bd0 (double x, 
     double np);

RESULT *
set_result_matrix(XDATA *xdata,
		  int *idx,
		  int hits);

void
calculate_pw_hyper(XDATA *xdata,
		   RESULT *res,
		   int hits,
		   char **ids,
		   double Q,
		   enum TAIL tail,
		   RESULT **resmatp,
		   RESULT **resmato,
		   int exist);

void
find_idx_for_ids(XDATA *xdata,
		 char **ids,
		 int hits, 
		 int *perm);

void
print_matrix_results(RESULT *resp,
		     RESULT *reso,
		     char *matrix_file);
#endif
