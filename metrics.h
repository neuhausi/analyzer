/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: metrics.h,v 1.9 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#ifndef metrics_h
#define metrics_h

#include <math.h>
#include <stdlib.h>
EXTERN int verbose;
EXTERN int debug;
#include "stats.h"
#include "utils.h"
#include "files.h" 

double
euclidean(VECTOR *v1, VECTOR *v2);

double
ieuclidean(VECTOR *v1, VECTOR *v2);

double
minkowski(VECTOR *v1, VECTOR *v2, double power);

double
iminkowski(VECTOR *v1, VECTOR *v2, double power);

double
manhattan(VECTOR *v1, VECTOR *v2);

double
imanhattan(VECTOR *v1, VECTOR *v2);

double
maximum(VECTOR *v1, VECTOR *v2);
 
double
imaximum(VECTOR *v1, VECTOR *v2);
 
double
pearson(VECTOR *v1, VECTOR *v2);
 
double
ipearson(VECTOR *v1, VECTOR *v2);

double
spearman(VECTOR *v1, VECTOR *v2);

double
ispearman(VECTOR *v1, VECTOR *v2);

double
lspearman(VECTOR *v1, VECTOR *v2);

double
ilspearman(VECTOR *v1, VECTOR *v2);

#endif
