/*********************************************************************
 Copyright 2006 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: function.h,v 1.8 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#ifndef function_h
#define function_h

#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <unistd.h>
#include <math.h>
EXTERN int verbose;
EXTERN int debug;
#include "utils.h"

void
fgauss(double x, double a[], double *y, double dyda[], int na);

void
nlls4p(double x, double a[], double *y, double dyda[], int na);

#endif

