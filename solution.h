/*********************************************************************
 Copyright 2004 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: solution.h,v 1.12 2009/03/20 21:58:27 neuhausi Exp $
**********************************************************************/

#ifndef solution_h
#define solution_h

#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <unistd.h>
#include <math.h>
EXTERN int verbose;
EXTERN int debug;
#include "utils.h"
#include "stats.h"
#include "function.h"

void
ludcmp(double **a, int n, int *indx, double *d);

void
lubksb(double **a, int n, int *indx, double *b);

void
svdcmp(double **a, int m, int n, double *w, double **v);

void
svbksb(double **u, double *w, double **v, int m, int n, double *b, double *x);

void
qrdcmp(double **a, int n, double *c, double *d, int *sing);

void
qrsolv(double **a, int n, double *c, double *d, double *b);

void
rsolv(double **a, int n, double *d, double *b);

void
qrupdt(double **r, double **qt, int n, double *u, double *v);

void
rotate(double **r, double **qt, int n, int i, double a, double b);

void
cholsl(double **a, int n, double *p, double *b, double *x);

void
choldc(double **a, int n, float *p);

void
cholroot(double **a, int n, float *p);

int
gaussj(double **a, int n, double **b, int m);

int
mrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[],
       int ma, double **covar, double **alpha, double *chisq,
       void (*funcs)(double, double [], double *, double [], int), double *alamda);

void
mrqcof(double x[], double y[], double sig[], int ndata, double a[], int ia[],
       int ma, double **alpha, double beta[], double **chisq,
       void (*funcs)(double, double [], double *, double [], int));

#endif

