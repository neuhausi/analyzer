/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

 $Id: methods.h,v 1.6 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#ifndef methods_h
#define methods_h

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
EXTERN int verbose;
EXTERN int debug;
#include "utils.h"
#include "files.h"

XDATA *
transpose_xdata(XDATA *xdata);

XDATA *
remove_null_xdata(XDATA *xdata,
		  FACTOR *factor,
		  int s);

MATRIX *
transpose_matrix(MATRIX *matrix);

MATRIX *
remove_null_matrix(MATRIX *matrix,
		   FACTOR *factor);

FACTOR *
remove_null_factor(FACTOR *factor,
		   int s);

VECTOR *
unique_factor(FACTOR *factor);

#endif
