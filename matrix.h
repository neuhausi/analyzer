/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

 $Id: matrix.h,v 1.16 2010/04/06 21:21:44 neuhausi Exp $
**********************************************************************/

#ifndef matrix_h
#define matrix_h

#include <stdlib.h>
#include <malloc.h>
#include <unistd.h>
#include <math.h>
EXTERN int verbose;
EXTERN int debug;
#include "utils.h"
#include "files.h"
#include "solution.h"

MATRIX *
mem_allocate_matrix(int n, int m);

MATRIX *
transpose(MATRIX *m,
	  int w);

double
dot(MATRIX *m1,
    MATRIX *m2,
    int mm1,
    int mm2);

MATRIX *
vector_to_matrix(VECTOR *v,
		 int w);

double
determinant(MATRIX *m,
	    int w);

MATRIX *
cofactor(MATRIX *m,
	 int t,
	 int w);

MATRIX * 
multiply1(MATRIX *m1,
	  MATRIX *m2,
	  int w);

MATRIX * 
multiply(MATRIX *m1,
	 MATRIX *m2,
	 int w);
	 
MATRIX * 
multiply_number(MATRIX *m1,
	  	double n,
	  	int w);

MATRIX *
diag_add(MATRIX *m,
	 VECTOR *v,
	 int w);

MATRIX * 
add(MATRIX *m1,
    MATRIX *m2,
    int w);

MATRIX * 
subtract(MATRIX *m1,
	 MATRIX *m2,
	 int w);

MATRIX *
inv(MATRIX *m,
    int w);

MATRIX *
diag(VECTOR *v);

MATRIX *
fill(int m,
     int n,
     double f);

MATRIX *
identity(int n);

MATRIX *
remove_cols(MATRIX *m,
	    int *filter,
	    int col);

MATRIX *
augment_cols(MATRIX *m1,
	     MATRIX *m2,
	     int w);

MATRIX *
remove_rows(MATRIX *m,
	    int *filter,
	    int row);

MATRIX *
augment_rows(MATRIX *m1,
	     MATRIX *m2,
	     int w);

MATRIX *
duplicate(MATRIX *m);

MATRIX *
copy(MATRIX *m,
     int n1,
     int m1,
     int n2,
     int m2);

MATRIX *
kronecker(MATRIX *m1,
	  MATRIX *m2,
          int w);

MATRIX *
cross(MATRIX *m1,
      MATRIX *m2,
      int w);

MATRIX *
clean_cross(MATRIX *m1,
	    int **zeros);

MATRIX *
invert(MATRIX *m,
       int w);

MATRIX *
pinvert(MATRIX *mat, VECTOR *vect);

MATRIX *
g2invert(MATRIX *m,
	 int k,
	 int **rank,
	 int *r,
	 double *det,
	 int w);
	 
void
sweep(MATRIX *m,
      int k1,
      int k2,
      int **rank,
      int *r,
      double *det);
      
void
sweep1(MATRIX *m,
      int k1,
      int k2,
      int **rank,
      int *r,
      double *det);  
          
void
mat_rank(MATRIX *m,
	 int **rank,
	 int *r,
	 double *det);

MATRIX *
independent(MATRIX *m);

double
trace_log(MATRIX *m,
	  int w);

double
trace(MATRIX *m,
      int w);

int
rank_matrix(MATRIX *m);

double
max_matrix(MATRIX *m);

double
min_matrix(MATRIX *m);

void
dump_matrix(MATRIX *mat);

MATRIX * 
extract_diag(MATRIX *m);

MATRIX *
matrix_power(MATRIX *m,
             double p,
             int w);

MATRIX *
matrix_element_manipulation(MATRIX *m1,
                            MATRIX *m2,
                            int type);
			         
#endif
