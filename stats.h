/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: stats.h,v 1.17 2007/10/08 16:24:41 neuhausi Exp $
**********************************************************************/

#ifndef stats_h
#define stats_h

#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <unistd.h>
#include <math.h>
EXTERN int verbose;
EXTERN int debug;
#include "utils.h"

#define MAXIT 100
#define MAX_ITER 2000
#define TOL 1.0e-7
#define EPS 1.0e-12
#define EPS1 1.0e-13
#define DBL_EPSILON 2.220446049250313e-16
#define EPSKS1 0.001
#define EPSKS2 1.0e-8
#define FPMIN 1.0e-30
#define TINY 1.0e-20
#define MD fmod(1,0)
#define PI 3.141593
#define PIX2 6.283185307179586476925286766559

#define SQR(a) ((a)*(a))

#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define MIN(a,b) ((a) > (b) ? (b) : (a))

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

#define SIGN(a, b) ((b) < 0 ? -fabs(a) : fabs(a))

void
simple_sort(int n, double *arr);

void
sort (int n, double *uar, double *sar);

void
sorti (int n, int *uar, int *sar);

void
sort_by_index (int n, double *uar, int *indx);

void
sort_by_indexi (int n, int *uar, int *indx);

void
sort_by_indexc (int n, char **uar, int *indx);

void
reverse(int n, int *arr);

void
pearsn(double x[], double y[], int n, double *r, double *prob,
       double *z);

void
spearman_corr(double data1[], double data2[], int n, int l, double *d, double *zd,
	      double *probd, double *rs, double *probrs);

void
spearman_corr_nna(double data1[], double data2[], int n, int l, double *d, double *zd,
		  double *probd, double *rs, double *probrs);

double
betai(double a, double b, double x);

double
betacf(double a, double b, double x);

double
gammln(double xx);

void
rank(int n, int indx[], int irank[]);

void
crank(int n, double w[], double *s);

void
crankl(int n, double w[], double *s);

double
erfcc(double x);

void
avevar(double data[], int n, double *ave, double *var);

void
bonferroni (double *data, double *cdata, int n);

void
fdr (int n, double *data, double *qval, double q,
     double eta0, double *cutoff, int *nsig);

void
histogram(double *data, int n, int intervals, double bin, double min, int *hist);

void
percentile(int *data, int n, int tot, double *prct);

void
shuffle(int *array, int n);

int
factrl (int n);

int
permutation (int n, int k);

int
choose (int n, int k);

double
pythag(double a, double b);

void
tqli(double d[], double e[], int n, double **z);

void
tred2(double **a, int n, double d[], double e[]);

void
covsrt(double **covar, int ma, int ia[], int mfit);

double
probks(double alam);

void
zscore(double **matrix,
       int n,
       int m,
       int axis);

#endif

