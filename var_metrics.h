/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net 

 $Id: var_metrics.h,v 1.8 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#ifndef var_metrics_h
#define var_metrics_h

#include <math.h>
#include <stdlib.h>
EXTERN int verbose;
EXTERN int debug;
#include "stats.h"
#include "utils.h"
#include "files.h" 

double
var_euclidean(VAR_VECTOR *v1, VECTOR *v2);

double
var_ieuclidean(VAR_VECTOR *v1, VECTOR *v2);

double
var_minkowski(VAR_VECTOR *v1, VECTOR *v2, double power);

double
var_iminkowski(VAR_VECTOR *v1, VECTOR *v2, double power);

double
var_manhattan(VAR_VECTOR *v1, VECTOR *v2);

double
var_imanhattan(VAR_VECTOR *v1, VECTOR *v2);

double
var_maximum(VAR_VECTOR *v1, VECTOR *v2);
 
double
var_imaximum(VAR_VECTOR *v1, VECTOR *v2);
 
double
var_pearson(VAR_VECTOR *v1, VECTOR *v2);
 
double
var_ipearson(VAR_VECTOR *v1, VECTOR *v2);

double
var_spearman(VAR_VECTOR *v1, VECTOR *v2);

double
var_ispearman(VAR_VECTOR *v1, VECTOR *v2);

double
var_lspearman(VAR_VECTOR *v1, VECTOR *v2);

double
var_ilspearman(VAR_VECTOR *v1, VECTOR *v2);

#endif
