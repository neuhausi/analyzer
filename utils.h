/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus and Robert Brucolleri
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

 $Id: utils.h,v 1.8 2007/04/13 22:20:27 neuhausi Exp $
**********************************************************************/

#ifndef utils_h
#define utils_h

#include <malloc.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
EXTERN int verbose;
EXTERN int debug;
#include "stats.h"

#define TYPE_ALLOC(type) (type *) mem_alloc(sizeof(type))

#define ARRAY_ALLOC(nelement, type) (type *) mem_alloc((nelement) * sizeof(type))

#define FREAD_CALL(call, n_elem, name) \
    if ((call) != (n_elem)) { \
        fprintf(stderr, "%s: Unable to read %d elements\n", name, n_elem); \
	exit(1); \
    } \

char
*mem_alloc(int n);

void
exit_on_error(char *str);

void
print_array_d(double data[], int n);

void
print_array_i(int data[], int n);

void
print_mat_d(double **data, int m, int n);

void
print_mat_i(int **data, int m, int n);

void
print_anova_data(double data[], int n, int **factor, int c);

#endif
