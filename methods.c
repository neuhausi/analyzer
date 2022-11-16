/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

 $Id: methods.c,v 1.10 2007/04/12 22:22:39 neuhausi Exp $
**********************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define EXTERN
#include "methods.h"

XDATA *
transpose_xdata(XDATA *xdata)
{

  int i, j, len;
  XDATA *x;
  VECTOR *v, *vx;
  
  x = TYPE_ALLOC(XDATA);
  x->n = xdata->s;
  x->s = xdata->n;
  x->md = ARRAY_ALLOC(x->n, int);
  x->variables = ARRAY_ALLOC(x->n, char *);
  x->samples = ARRAY_ALLOC(x->s, char *);
  x->vectors = ARRAY_ALLOC(x->n, VECTOR *);
  
  for (i = 0; i < x->s; i++) {
    len = strlen(xdata->variables[i]);
    x->samples[i] = ARRAY_ALLOC((len+1), char);
    x->samples[i] = xdata->variables[i];
    x->samples[i][len] = '\0';
  }
  
  for (i = 0; i < x->n; i++) {
    len = strlen(xdata->samples[i]);
    x->variables[i] = ARRAY_ALLOC((len+1), char);
    x->variables[i] = xdata->samples[i];
    x->variables[i][len] = '\0';
  }
  
  for (i = 0; i < x->n; i++) {
    x->md[i] = 0;
    vx = x->vectors[i] = TYPE_ALLOC(VECTOR);
    vx->data = ARRAY_ALLOC(x->s, double);
    vx->n = x->s;
    for (j = 0; j < x->s; j++) {
      v = xdata->vectors[j];
      vx->data[j] = v->data[i];
      if (isnan(vx->data[j])) {
	x->md[i]++;  
      }
    }
  }
  
  free(xdata);
  return x;
  
}

XDATA *
remove_null_xdata(XDATA *xdata,
		  FACTOR *factor,
		  int s)
{
 
  // s is the new sample length
 
  int i, j, k, len;
  XDATA *x;
  VECTOR *v, *vn;

  x = TYPE_ALLOC(XDATA);
  x->n = xdata->n;
  x->s = s;
  x->md = ARRAY_ALLOC(x->n, int);
  x->variables = ARRAY_ALLOC(x->n, char *);
  x->samples = ARRAY_ALLOC(x->s, char *);
  x->vectors = ARRAY_ALLOC(x->n, VECTOR *);
  
  // The variables names
  for (i = 0; i < xdata->n; i++) {
    len = strlen(xdata->variables[i]);
    x->variables[i] = ARRAY_ALLOC((len+1), char);
    x->variables[i] = xdata->variables[i];
    x->variables[i][len] = '\0';
  }
  
  // The sample names
  k = 0;
  for (i = 0; i < xdata->s; i++) {
    if (factor->factors[0][i] >= 0) {
      len = strlen(xdata->samples[i]);
      x->samples[k] = ARRAY_ALLOC((len+1), char);
      x->samples[k] = xdata->samples[i];
      x->samples[k][len] = '\0';
      k++;
    }
  }
  
  // The data
  for (i = 0; i < xdata->n; i++) {
    
    x->md[i] = 0;

    vn = xdata->vectors[i];
    v = x->vectors[i] = TYPE_ALLOC(VECTOR);
    v->data = ARRAY_ALLOC(s, double);
    v->n = s;
   
    k = 0;
    for (j = 0; j < xdata->s; j++) {
      if (factor->factors[0][i] >= 0) {
	x->vectors[i]->data[k] = vn->data[j];
	if (isnan(x->vectors[i]->data[k])) {
	  x->md[i]++;  
	}
	k++;
      }
    }
    
  }  
  free(xdata);

  return x;
  
}

FACTOR *
remove_null_factor(FACTOR *factor,
		   int s)
{

  // s is the new sample length

  int i, j, k;
  
  FACTOR *f;
  
  f = TYPE_ALLOC(FACTOR);
  f->n = factor->n;
  f->s = s;
  f->classes = ARRAY_ALLOC(f->n, int);
  f->members = ARRAY_ALLOC(f->n, int *);
  f->factors = ARRAY_ALLOC(f->n, int *);
  
  for (i = 0; i < f->n; i++) {
    
    f->factors[i] = ARRAY_ALLOC(s, int);
    f->members[i] = ARRAY_ALLOC(factor->classes[i], int);
    f->classes[i] = factor->classes[i];
    
    for (j = 0; j < factor->classes[i]; j++) {
      f->members[i][j] = factor->members[i][j];
    }
    
    k = 0;
    for (j = 0; j < factor->s; j++) {
      if (factor->factors[0][j] >= 0) {
	factor->factors[i][k] = factor->factors[i][j];
	k++;
      }
    }
  }
  free(factor);

  return f;
  
}

VECTOR *
unique_factor(FACTOR *factor)
{
  
  int i, j, k, cnt, mult;
  int uniques[1000000];
  VECTOR *v;
  
  v = TYPE_ALLOC(VECTOR);
  v->data = ARRAY_ALLOC(factor->n, double);
  v->n = factor->n;
  
  for (i = 0; i < pow(10, (factor->n+1)); i++) {
    uniques[i] = -1;
  }
  
  cnt = 0;
  for (i = 0; i < factor->s; i++) {
    mult = 1;
    for (j = 0; j < factor->n; j++) {
      k += (factor->factors[j][i] * mult);
      mult *= 10;
    }
    if (uniques[k] == -1) {
      uniques[k] = cnt;
      v->data[i] = uniques[k];
      cnt++;
    } else {
      v->data[i] = uniques[k];
    }
  }
  
  return v;

}

MATRIX *
transpose_matrix(MATRIX *matrix)
{
  
  int i, j, len;
  MATRIX *x;
  
  x = TYPE_ALLOC(MATRIX);
  x->n = matrix->m;
  x->m = matrix->n;
  x->md = ARRAY_ALLOC(x->n, int);
  x->rows = ARRAY_ALLOC(x->n, char *);
  x->cols = ARRAY_ALLOC(x->m, char *);
  x->data = ARRAY_ALLOC(x->n, double *);
  
  for (i = 0; i < x->m; i++) {
    len = strlen(matrix->rows[i]);
    x->cols[i] = ARRAY_ALLOC((len+1), char);
    x->cols[i] = matrix->rows[i];
    x->cols[i][len] = '\0';
  }
  
  for (i = 0; i < x->n; i++) {
    len = strlen(matrix->cols[i]);
    x->rows[i] = ARRAY_ALLOC((len+1), char);
    x->rows[i] = matrix->cols[i];
    x->rows[i][len] = '\0';
  }
  
  for (i = 0; i < x->n; i++) {
    x->md[i] = 0;
    x->data[i] = ARRAY_ALLOC(x->m, double);
    for (j = 0; j < x->m; j++) {
      x->data[i][j] = matrix->data[j][i];
      if (isnan(x->data[i][j])) {
	x->md[i]++;  
      }
    }
  }  
  free(matrix);

  return x;
  
}

MATRIX *
remove_null_matrix(MATRIX *matrix,
		   FACTOR *factor)
{
  
  // s is the new sample length

  int i, j, k, len;
  MATRIX *x;
  
  x = TYPE_ALLOC(MATRIX);
  x->n = matrix->n;
  x->md = ARRAY_ALLOC(x->n, int);
  x->rows = ARRAY_ALLOC(x->n, char *);
  x->cols = ARRAY_ALLOC(x->m, char *);
  x->data = ARRAY_ALLOC(x->n, double *);
  
  // The variables names
  for (i = 0; i < matrix->n; i++) {
    len = strlen(matrix->rows[i]);
    x->rows[i] = ARRAY_ALLOC((len+1), char);
    x->rows[i] = matrix->rows[i];
    x->rows[i][len] = '\0';
  }
  
  // The sample names
  k = 0;
  for (i = 0; i < matrix->m; i++) {
    if (factor->factors[0][i] >= 0) {
      len = strlen(matrix->cols[i]);
      x->cols[k] = ARRAY_ALLOC((len+1), char);
      x->cols[k] = matrix->cols[i];
      x->cols[k][len] = '\0';
      k++;
    }
  }
  x->m = k;
  
  // The data
  for (i = 0; i < matrix->n; i++) {
    x->md[i] = 0;
    x->data[i] = ARRAY_ALLOC(x->m, double);
    k = 0;
    for (j = 0; j < matrix->m; j++) {
      if (factor->factors[0][i] >= 0) {
	x->data[i][k] = matrix->data[i][j];
	if (isnan(x->data[i][k])) {
	  x->md[i]++;  
	}
	k++;
      }
    }
  }
  free(matrix);

  return x;
  
}
