/*********************************************************************
 Copyright 2006 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: function.c,v 1.9 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <unistd.h>
#include <math.h>
#define EXTERN
#include "function.h"

void
fgauss(double x, double a[], double *y, double dyda[], int na)
{
  int i;
  double fac,ex,arg;

  *y=0.0;
  for (i=0;i<na-1;i+=3) {
    arg=(x-a[i+1])/a[i+2];
    ex=exp(-arg*arg);
    fac=a[i]*ex*2.0*arg;
    *y += a[i]*ex;
    dyda[i]=ex;
    dyda[i+1]=fac/a[i+2];
    dyda[i+2]=fac*arg/a[i+2];
  }

}

void
nlls4p(double x, double a[], double *y, double dyda[], int na)
{

  if (na != 4) {
    exit_on_error("not the expected number of parameters.");
  }

  double dose = pow(x, a[3]);

  *y = a[1] + (a[0] - a[1])/(1.0 + a[2]/dose);

  dyda[0] = 1.0 / (1.0 + a[2] / dose);
  dyda[1] = 1.0 - dyda[0];
  dyda[2] = -(a[0] - a[1]) * dyda[0] * dyda[0] / dose;
  dyda[3] = (a[0] - a[1]) * dyda[0] * dyda[0] * log(x) * a[2] / dose;

}
