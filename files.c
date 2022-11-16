/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: files.c,v 1.34 2008/07/02 22:16:14 neuhausi Exp $
**********************************************************************/

#include <malloc.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#define EXTERN
#include "files.h"

int
is_binary_file(char *file_name)
{

  /* This subrutine was contributed by John Hinsdale */

  FILE *in;
  int is_binary = 0;
  int maxread = 10000;
  int nread = 0;
  int a_char = 0;

  if ((in = fopen(file_name, "r")) == NULL) { 
    fprintf(stderr, "Error in opening file %s\n", file_name);
    perror(0);
    exit(1);
  }


  while ( EOF != (a_char = fgetc(in)) ) {
    nread++;
    if (nread >= maxread) {
      break;
    }
    /* Allow TAB, CR, LF, FF, etc. */
    if (a_char < 0x09 ||
	(a_char > 0x0d && a_char < 0x20)
	|| a_char > 0xff) {
      is_binary = 1;
      break;
    }
  }

  fclose(in);

  return is_binary;

}

int
file_exists(char *file_name)
{

  FILE *in;
  int exist;

  if ((in = fopen(file_name, "r")) == NULL) { 
    exist = 0;
  } else {
    fclose(in);
    exist = 1;
  }

  return exist;

}

XDATA *
read_xdata(char *file_name, 
	   int n_samples,
	   int n_variables)
{
  FILE *in;
  int i, j, tlen;
  VECTOR *v;
  XDATA *x;
  char line[8000000];
  char *tok;
  
  x = TYPE_ALLOC(XDATA);
  x->n = n_variables;
  x->s = n_samples;
  x->md = ARRAY_ALLOC(x->n, int);
  x->variables = ARRAY_ALLOC(x->n, char *);
  x->samples = ARRAY_ALLOC(x->s, char *);
  x->vectors = ARRAY_ALLOC(x->n, VECTOR *);
  
  if ((in = fopen(file_name, "r")) == NULL) { 
    fprintf(stderr, "Error in opening file %s\n", file_name);
    perror(0);
    exit(1);
  }
  
  /* Read the header */
  
  if (fgets(line, sizeof(line), in) == NULL) {
    fprintf(stderr, "Premature end of file %s\n", file_name);
    exit(1);
  }
  
  i = 0;
  tok = strtok (line, "\t\n");
  while (tok != NULL) {
    tlen = strlen(tok);
    x->samples[i] = ARRAY_ALLOC((tlen+1), char);
    strncpy(x->samples[i], tok, tlen);
    x->samples[i][tlen] = '\0';
    i++;
    tok = strtok (NULL, "\t\n");
  }
  if (i != n_samples) {
    fprintf(stderr, "Bad sample name count in file %s (%i / %i)\n", file_name, i, n_samples);
    exit(1);
  }
  
  /* Read the rest of the file */
  
  for (i = 0; i < n_variables; i++) {
    
    /* No missing data */
    x->md[i] = 0;

    v = x->vectors[i] = TYPE_ALLOC(VECTOR);
    v->data = ARRAY_ALLOC(n_samples, double);
    v->n = n_samples;
    
    if (fgets(line, sizeof(line), in) == NULL) {
      fprintf(stderr, "Premature end of file %s\n", file_name);
      exit(1);
    }
    
    tok = strtok (line, "\t\n");
    tlen = strlen(tok);
    x->variables[i] = ARRAY_ALLOC((tlen+1), char);
    strncpy(x->variables[i], tok, tlen);
    x->variables[i][tlen] = '\0';
    tok = strtok (NULL, "\t\n");
    j = 0;
    
    while (tok != NULL) {
      if (!strcmp(tok, "NA") == 0) {
	v->data[j] = (double)atof(tok);
      } else {
	v->data[j] = MD;
	x->md[i]++;
      }
      j++;
      tok = strtok (NULL, "\t\n");
    }
    if (j != n_samples) {
	fprintf(stderr, "Bad sample name count in file %s (%i / %i) in line %i %s\n", file_name, j, n_samples, i, x->samples[i]);
      exit(1);
    }

  }
  
  fclose(in);

  return x;
    
}

XDATA *
read_binary_xdata(char *file_name)
{
  FILE *in;
  int n_variables, n_samples, i, j, len;
  VECTOR *v;
  XDATA *x;
  
  if ((in = fopen(file_name, "r")) == NULL) { 
    fprintf(stderr, "Error in opening file %s\n", file_name);
    perror(0);
    exit(1);
  }
  
  FREAD_CALL(fread(&n_variables, sizeof(n_variables), 1, in), 1, "n_variables");
  FREAD_CALL(fread(&n_samples, sizeof(n_samples), 1, in), 1, "n_samples");
  
  x = TYPE_ALLOC(XDATA);
  x->n = n_variables;
  x->s = n_samples;
  x->md = ARRAY_ALLOC(x->n, int);
  x->variables = ARRAY_ALLOC(x->n, char *);
  x->samples = ARRAY_ALLOC(x->s, char *);
  x->vectors = ARRAY_ALLOC(x->n, VECTOR *);
  
  for (i = 0; i < x->s; i++) {
    FREAD_CALL(fread(&len, sizeof(len), 1, in), 1, "sample name length");
    x->samples[i] = ARRAY_ALLOC((len + 1), char);
    FREAD_CALL(fread(x->samples[i], len, 1, in), 1, "sample name");
    x->samples[i][len] = '\0';
  }
  
  for (i = 0; i < x->n; i++) {
    
    x->md[i] = 0;

    v = x->vectors[i] = TYPE_ALLOC(VECTOR);
    v->data = ARRAY_ALLOC(x->s, double);
    v->n = x->s;
    
    FREAD_CALL(fread(&len, sizeof(len), 1, in), 1, "variable name length");
    x->variables[i] = ARRAY_ALLOC((len + 1), char);
    FREAD_CALL(fread(x->variables[i], len, 1, in), 1, "variable name");
    x->variables[i][len] = '\0';
    FREAD_CALL(fread(v->data, sizeof(double), x->s, in), x->s, "data vector");

    for (j = 0; j < x->s; j++) {
      if (isnan(v->data[j])) {
	x->md[i]++;  
      }  
    }

  }
  
  fclose(in);

  return x;
  
}

MATRIX *
read_matrix(char *file_name, 
	     int m,
	     int n)
{
  FILE *in;
  int i, j, tlen;
  MATRIX *x;
  char line[8000000];
  char *tok;
  
  x = TYPE_ALLOC(MATRIX);
  x->n = n;
  x->m = m;
  x->md = ARRAY_ALLOC(x->n, int);
  x->rows = ARRAY_ALLOC(x->n, char *);
  x->cols = ARRAY_ALLOC(x->m, char *);
  x->data = ARRAY_ALLOC(x->n, double *);
  
  if ((in = fopen(file_name, "r")) == NULL) { 
    fprintf(stderr, "Error in opening file %s\n", file_name);
    perror(0);
    exit(1);
  }
  
  /* Read the header */
  
  if (fgets(line, sizeof(line), in) == NULL) {
    fprintf(stderr, "Premature end of file %s\n", file_name);
    exit(1);
  }
  
  i = 0;
  tok = strtok (line, "\t\n");
  while (tok != NULL) {
    tlen = strlen(tok);
    x->cols[i] = ARRAY_ALLOC((tlen+1), char);
    strncpy(x->cols[i], tok, tlen);
    x->cols[i][tlen] = '\0';
    i++;
    tok = strtok (NULL, "\t\n");
  }
  if (i != m) {
    fprintf(stderr, "Bad sample name count in file %s\n", file_name);
    exit(1);
  }
  
  /* Read the rest of the file */
  
  for (i = 0; i < n; i++) {
    
    /* No missing data */
    x->md[i] = 0;

    x->data[i] = ARRAY_ALLOC(m, double);
    
    if (fgets(line, sizeof(line), in) == NULL) {
      fprintf(stderr, "Premature end of file %s\n", file_name);
      exit(1);
    }
    
    tok = strtok (line, "\t\n");
    tlen = strlen(tok);
    x->rows[i] = ARRAY_ALLOC((tlen+1), char);
    strncpy(x->rows[i], tok, tlen);
    x->rows[i][tlen] = '\0';
    tok = strtok (NULL, "\t\n");
    j = 0;
    
    while (tok != NULL) {
      if (!strcmp(tok, "NA") == 0) {
	x->data[i][j] = (double)atof(tok);
      } else {
	x->data[i][j] = MD;
	x->md[i]++;
      }
      j++;
      tok = strtok (NULL, "\t\n");
    }
    if (j != m) {
      fprintf(stderr, "Bad sample name count in file %s\n", file_name);
      exit(1);
    }
  }
  
  fclose(in);

  return x;
  
}

MATRIX *
read_binary_matrix(char *file_name)
{
  FILE *in;
  int n, m, i, j, len;
  char *name;
  MATRIX *x;
  
  if ((in = fopen(file_name, "r")) == NULL) { 
    fprintf(stderr, "Error in opening file %s\n", file_name);
    perror(0);
    exit(1);
  }

  FREAD_CALL(fread(&n, sizeof(n), 1, in), 1, "n");
  FREAD_CALL(fread(&m, sizeof(m), 1, in), 1, "m");

  x = TYPE_ALLOC(MATRIX);
  x->n = n;
  x->m = m;
  x->md = ARRAY_ALLOC(x->n, int);
  x->rows = ARRAY_ALLOC(x->n, char *);
  x->cols = ARRAY_ALLOC(x->m, char *);
  x->data = ARRAY_ALLOC(x->n, double *);
 
  for (i = 0; i < x->m; i++) {
    FREAD_CALL(fread(&len, sizeof(len), 1, in), 1, "sample name length");
    name = ARRAY_ALLOC((len + 1), char);
    FREAD_CALL(fread(name, len, 1, in), 1, "sample name");
    name[len] = '\0';
    x->cols[i] = name;
  }
  
  for (i = 0; i < x->n; i++) {	

    x->md[i] = 0;

    x->data[i] = ARRAY_ALLOC(x->m, double);

    FREAD_CALL(fread(&len, sizeof(len), 1, in), 1, "variable name length");
    name = ARRAY_ALLOC((len + 1), char);
    FREAD_CALL(fread(name, len, 1, in), 1, "variable name");
    name[len] = '\0';
    x->rows[i] = name;
    FREAD_CALL(fread(x->data[i], sizeof(double), x->m, in), x->m, "data matrix");

    for (j = 0; j < x->m; j++) {
      if (isnan(x->data[i][j])) {
	x->md[i]++;  
      }  
    }

  }

  fclose(in);

  free(name);

  return x;
    
}

VECTOR *
read_pattern(char *file_name, 
	     int n)
{
  
  FILE *in;
  int i;
  VECTOR *v;
  char line[8000000];
  char *tok;
  
  v = TYPE_ALLOC(VECTOR);
  v->data = ARRAY_ALLOC(n, double);
  v->n = n;
  
  if ((in = fopen(file_name, "r")) == NULL) { 
    fprintf(stderr, "Error in opening file %s\n", file_name);
    perror(0);
    exit(1);
  }
  
  // This line will contain the sample names but we don't check it
  if (fgets(line, sizeof(line), in) == NULL) {
    fprintf(stderr, "Premature end of file %s\n", file_name);
    exit(1);
  }

  if (fgets(line, sizeof(line), in) == NULL) {
    fprintf(stderr, "Premature end of file %s\n", file_name);
    exit(1);
  }
  
  i = 0;
  tok = strtok (line, "\t");
  while (tok != NULL) {
    if (!strcmp(tok, "NA") == 0) {
      v->data[i] = (double)atof(tok);
    } else {
      v->data[i] = MD;
    }
    i++;
    tok = strtok (NULL, "\t\n");
  }
  if (i != n) {
    fprintf(stderr, "Bad sample name count in file %s\n", file_name);
    exit(1);
  }
  
  fclose(in);
  
  return v;
  
}

FACTOR *
read_factor(char *file_name, 
	    int n_factors,
            int n_samples)
{
  
  FILE *in;
  int i, j, k, seen, cnt, tlen;
  int *c, *m;
  int array[n_samples]; 
  int classes[n_samples]; 
  int n_class[n_samples];
  int members[n_samples];
  FACTOR *f;
  
  char line[8000000];
  char *tok, **n;
  char *temp[10000];
  char *names[10000];
  
  f = TYPE_ALLOC(FACTOR);
  f->n = n_factors;
  f->s = n_samples;

  if (n_factors==0) {
    return f;
  }

  f->classes = ARRAY_ALLOC(f->n, int);
  f->members = ARRAY_ALLOC(f->n, int *);
  f->factors = ARRAY_ALLOC(f->n, int *);
  f->samples = ARRAY_ALLOC(f->s, char *);
  f->fnames  = ARRAY_ALLOC(f->n, char *);
  f->cnames  = ARRAY_ALLOC(f->n, char **);
 
  if ((in = fopen(file_name, "r")) == NULL) { 
    fprintf(stderr, "Error in opening file %s\n", file_name);
    perror(0);
    exit(1);
  }
  
  /* Read the header */
  
  if (fgets(line, sizeof(line), in) == NULL) {
    fprintf(stderr, "Premature end of file %s\n", file_name);
    exit(1);
  }

  i = 0;
  tok = strtok (line, "\t\n");
  while (tok != NULL) {
    tlen = strlen(tok);
    f->samples[i] = ARRAY_ALLOC((tlen+1), char);
    strncpy(f->samples[i], tok, tlen);
    f->samples[i][tlen] = '\0';
    i++;
    tok = strtok (NULL, "\t\n");
  }

  if (i != n_samples) {
    fprintf(stderr, "Bad sample name count in file %s\n", file_name);
    exit(1);
  }
  
  /* Read the rest of the file */
  
  for (i = 0; i < n_factors; i++) {
    
    if (fgets(line, sizeof(line), in) == NULL) {
      fprintf(stderr, "Premature end of file %s\n", file_name);
      exit(1);
    }

    /* Get the factor name */    
    
    tok = strtok (line, "\t\n");
    tlen = strlen(tok);
    f->fnames[i] = ARRAY_ALLOC((tlen+1), char);
    strncpy(f->fnames[i], tok, tlen);
    f->fnames[i][tlen] = '\0';
    tok = strtok (NULL, "\t\n");

    /* Parse the rest of the line */
    /* and put in temporary array */
        
   /* Reset values for all arrays */
        
    for (j = 0; j < n_samples; j++) {
      array[j]   = 0; 
      classes[j] = 0; 
      n_class[j] = 0;
      members[j] = 0;
    }

    j = 0;
    cnt = 0;
    while (tok != NULL) {
      tlen = strlen(tok);
      temp[j] = ARRAY_ALLOC((tlen+1), char);
      strncpy(temp[j], tok, tlen);
      temp[j][tlen] = '\0';
      if (strcmp(temp[j], "Unassigned") == 0) {
	array[j] = -1;   
      } else {
	if (j > 0) {
	  seen = 0;
	  for (k = 0; k < j; k++) {
	    if (strcmp(temp[j], temp[k]) == 0){
	      array[j] = array[k];
	      seen = 1;
	    }	    
	  }
	  if (seen == 0) {
	    array[j] = cnt;
	    cnt++;
	  }
	} else {
	  array[j] = cnt;
	  cnt++;
	}
      }
      j++;
      tok = strtok (NULL, "\t\n");
    }

    if (j != n_samples) {
      fprintf(stderr, "Bad sample name count in file %s\n", file_name);
      exit(1);
    }
        
    /* Assign the classes in each factor and make
       sure that the classes go from 0 to n - 1 */
    
    c = f->factors[i] = ARRAY_ALLOC(n_samples, int);

    cnt = 0;

    for (j = 0; j < n_samples; j++) {
      if (array[j] >= 0 ) {
	if (classes[array[j]] == 0) {
	  classes[array[j]]++;
	  n_class[array[j]] = cnt;
	  members[cnt]++;
	  c[j] = cnt;
	  tlen = strlen(temp[j]);
	  names[cnt] = ARRAY_ALLOC((tlen+1), char);
	  strncpy(names[cnt], temp[j], tlen);
	  names[cnt][tlen] = '\0';
	  cnt++;
	} else {
	  members[n_class[array[j]]]++;
	  c[j] = n_class[array[j]];
	}
      } else {
	c[j] = array[j];
      }
    }

    f->classes[i] = cnt;
    m = f->members[i] = ARRAY_ALLOC(cnt, int);
    n = f->cnames[i] = ARRAY_ALLOC(cnt, char *);
   
    for (j = 0; j < cnt; j++) {
      m[j] = members[j];
      tlen = strlen(names[j]);
      n[j] = ARRAY_ALLOC((tlen+1), char);
      strncpy(n[j], names[j], tlen);
      n[j][tlen] = '\0';
    }

  }
  
  fclose(in);

  if (debug) {
    printf ("FACTOR\n");
    printf ("\n");
    printf (" Number of Samples: %i\n", f->s);
    printf ("\n");
    printf (" Number of Factors: %i\n", f->n);
    printf ("\n");
    printf (" Factors\n");
    printf ("   S  ");
    for (i = 0; i < f->n; i++) {
      printf (" F%-3i", i);
    }
    printf ("\n");
    for (i = 0; i < f->s; i++) {    
      printf (" %3i ", i+1);
      for (j = 0; j < f->n; j++) {
	printf ("%4i ", f->factors[j][i]);
      }
      printf ("\n");
    }
    printf ("\n");
    printf (" Classes in each factor\n");
    printf (" ");
    print_array_i(f->classes, f->n);
    printf ("\n");
    printf (" Members in each class:\n");
    for (i = 0; i < f->n; i++) {
      printf (" F%-2i=> ", i);
      print_array_i(f->members[i], f->classes[i]);
    }
    printf ("\n");
  }   

  return f;

}

COVARIATE *
read_covariate(char *file_name, 
	       int n_covariates,
	       int n_samples)
{
  
  FILE *in;
  int i, j, tlen;
  COVARIATE *c;

  char line[8000000];
  char *tok;  
  char *temp;

  c = TYPE_ALLOC(COVARIATE);
  c->n = n_covariates;

  if (n_covariates==0) {
    return c;
  }

  c->s = n_samples;
  c->samples = ARRAY_ALLOC(c->s, char *);
  c->cnames  = ARRAY_ALLOC(c->n, char *);
  c->values  = ARRAY_ALLOC(c->n, double *);
 
  if ((in = fopen(file_name, "r")) == NULL) { 
    fprintf(stderr, "Error in opening file %s\n", file_name);
    perror(0);
    exit(1);
  }
  
  /* Read the header */

  if (fgets(line, sizeof(line), in) == NULL) {
    fprintf(stderr, "Premature end of file %s\n", file_name);
    exit(1);
  }
  
  i = 0;
  tok = strtok (line, "\t\n");
  while (tok != NULL) {
    tlen = strlen(tok);
    c->samples[i] = ARRAY_ALLOC((tlen+1), char);
    strncpy(c->samples[i], tok, tlen);
    c->samples[i][tlen] = '\0';
    i++;
    tok = strtok (NULL, "\t\n");
  }

  if (i != n_samples) {
    fprintf(stderr, "Bad sample name count in file %s\n", file_name);
    exit(1);
  }

  /* Read the rest of the file */
  
  for (i = 0; i < n_covariates; i++) {

    if (fgets(line, sizeof(line), in) == NULL) {
      fprintf(stderr, "Premature end of file %s\n", file_name);
      exit(1);
    }

    /* Get the covariate name */    
    
    tok = strtok (line, "\t\n");
    tlen = strlen(tok);
    c->cnames[i] = ARRAY_ALLOC((tlen+1), char);
    strncpy(c->cnames[i], tok, tlen);
    c->cnames[i][tlen] = '\0';
    tok = strtok (NULL, "\t\n");

    /* Parse the rest of the line */

    c->values[i] = ARRAY_ALLOC(c->s, double);

    j = 0;
    while (tok != NULL) {
      tlen = strlen(tok);
      temp = ARRAY_ALLOC((tlen+1), char);
      strncpy(temp, tok, tlen);
      temp[tlen] = '\0';
      if (strcmp(temp, "Unassigned") == 0) {
	c->values[i][j] = MD;
      } else {
	c->values[i][j] = (double)atof(tok);
      }
      j++;
      tok = strtok (NULL, "\t\n");
    }

    if (j != n_samples) {
      fprintf(stderr, "Bad sample name count in file %s\n", file_name);
      exit(1);
    }

  }
  
  fclose(in);
  
  free(temp);

  return c;
  
}

GROUP *
set_group (FACTOR *factor,
	   COVARIATE *covariate,
	   char *fit_factor,
	   char *agg_factor,
	   char *zero,
	   int logtrans)
{

  int i, ii, j, z, p, m, len;
  char str[1000];

  GROUP *g;
  
  z = -1;
  g = TYPE_ALLOC(GROUP);
 
  // Find the covariate number

  g->c = -1;
  if (covariate->n > 0) {
    for (i = 0; i < covariate->n; i++) {
      if (strcmp(covariate->cnames[i], fit_factor) == 0) {
	g->c = i;
      }
    }
  }
  if (g->c < 0) {
    fprintf(stderr, "Unable to find %s in covariate file\n", fit_factor);
    exit(1);
  }
  
  // We create one group if there is no factor object or create
  // groups for each one of the classes in the selected factor

  g->n = -1;
  for (i = 0; i < factor->n; i++) {
    if (strcmp(factor->fnames[i], agg_factor) == 0) {
      // Find the factor number
      g->f = i;
      g->n = factor->classes[i];
      break;
    }
  }
  if (g->n < 0) {
    fprintf(stderr, "Unable to find %s in factor file\n", agg_factor);
    exit(1);
  }

  // We find the level to be used as zero if it is defined
  if (zero != NULL) {
    for (i = 0; i < factor->classes[g->f]; i++) {
      // Find the level number
      if (strcmp(factor->cnames[g->f][i], zero) == 0) {
	z = i;
      }	
    }    
    if (z < 0) {
      fprintf(stderr, "Unable to find %s among levels for factor %s\n", zero, agg_factor);
      exit(1);
    }
    g->n--;
  }

  // Allocate memory to store sample indeces for the number of groups
  g->samples = ARRAY_ALLOC(g->n, int);
  g->levels = ARRAY_ALLOC(g->n, int);
  g->groups = ARRAY_ALLOC(g->n, int *);
  g->names = ARRAY_ALLOC(g->n, char *);

  // Allocate memory to store the y values for the covariate
  g->values = ARRAY_ALLOC(g->n, double *);

  if (z < 0) {

    // There is no zero level
    for (i = 0; i < g->n; i++) {

      p = factor->n > 0 ? factor->members[g->f][i] : covariate->s;
      g->samples[i] = p;

      // Allocate memory for the number of members in each class 
      g->groups[i] = ARRAY_ALLOC(p, int);

      // Allocate memory for the covariate value of each member in each class 
      g->values[i] = ARRAY_ALLOC(p, double);

      // Set the name of the level
      if (factor->n > 0) {
	strcpy(str, factor->cnames[g->f][i]);
      } else {
	strcpy(str, "fit");
      }
      len = strlen(str);
      g->names[i] = ARRAY_ALLOC((len+1), char);  
      strncpy(g->names[i], str, len);
      g->names[i][len] = '\0';
 
      // Set the data
      m=0;
      for (j = 0; j < covariate->s; j++) {
	if (i == factor->factors[g->f][j]) {
	  g->groups[i][m] = j;
	  g->values[i][m] = covariate->values[g->c][j];
	  if (logtrans > 0) {
	    g->values[i][m] = log(g->values[i][m]);
	  }
	  m++;
	  g->levels[i] = m;
	}
      }

    }

  } else {

    ii = 0;

    for (i = 0; i < g->n + 1; i++) {

      if (i != z) {

        p = factor->n > 0 ? factor->members[g->f][i] + factor->members[g->f][z] : covariate->s;
        g->samples[ii] = p;
        
        // Allocate memory for the number of members in each class 
        g->groups[ii] = ARRAY_ALLOC(p, int);
        
        // Allocate memory for the covariate value of each member in each class 
        g->values[ii] = ARRAY_ALLOC(p, double);
        
        // Set the name of the level
        if (factor->n > 0) {
          strcpy(str, factor->cnames[g->f][i]);
        } else {
          strcpy(str, "fit");
        }
        len = strlen(str);
        g->names[ii] = ARRAY_ALLOC((len+1), char);  
        strncpy(g->names[ii], str, len);
        g->names[ii][len] = '\0';
        
        // Set the data
        m=0;
        for (j = 0; j < covariate->s; j++) {
          if (i == factor->factors[g->f][j] || z == factor->factors[g->f][j]) {
            g->groups[ii][m] = j;
            g->values[ii][m] = covariate->values[g->c][j];
            if (logtrans > 0 && g->values[ii][m] > 0) {
              g->values[ii][m] = log(g->values[ii][m]);
            }
            m++;
            g->levels[ii] = m;
          }
        }

	ii++;

      }

    }

  }

  return g;

}

char **
read_file_column(char *file_name,
		 int rows,
                 int header)
{

  FILE *in;
  int i, j, tlen;
  char line[8000000];
  char *tok, **objects;

  objects = ARRAY_ALLOC(rows, char *);

  if ((in = fopen(file_name, "r")) == NULL) { 
    fprintf(stderr, "Error in opening file %s\n", file_name);
    perror(0);
    exit(1);
  }

  if (header) {
      
    /* We skip the line */
    if (fgets(line, sizeof(line), in) == NULL) {
      fprintf(stderr, "Premature end of file %s\n", file_name);
      exit(1);
    }

  }

  /* Read the rest of the file */

  for (i = 0; i < rows; i++) {

    if (fgets(line, sizeof(line), in) == NULL) {
      fprintf(stderr, "Premature end of file %s\n", file_name);
      exit(1);
    }

    tok = strtok (line, "\t\n");
    tlen = strlen(tok);
    objects[i] = ARRAY_ALLOC((tlen+1), char);
    strncpy(objects[i], tok, tlen);
    objects[i][tlen] = '\0';
    tok = strtok (NULL, "\t\n");
    j = 0;

  }

  return objects;

}

void
print_results(char *out_file,
	      RESULT *res,
	      int order,
	      int indx)
{
  
  // order
  // 0 : no order
  // 1 : ascending
  // 2 : descending

  FILE *out;
  char *name;
  int i, j, len, idx, perm[res->n];
  double tmp[res->n];
  VECTOR *v;

  if (res->n > 1) {
 
   if (order != 0) {
      for (i = 0; i < res->n; i++) {
	v = res->vectors[i];
	tmp[i] = v->data[indx];		
      }
      sort_by_index(res->n, tmp, perm);
      if (order == 2) {
	reverse(res->n, perm);
      }
    } else {
      for (i = 0; i < res->n; i++) {
	perm[i] = i;
      }
    }

  } else {

    perm[0]=0;

  }
    
  if (out_file == NULL) {
    
    for (i = 0; i < res->a; i++) {
      printf("\t%s", res->attributes[i]);
    }
    printf("\n");
    
    for (i = 0; i < res->n; i++) {
      
      idx = perm[i];
      printf("%s", res->variables[idx]);
      
      v = res->vectors[idx];
      for (j = 0; j < res->a; j++) {
	printf("\t%e", v->data[j]);
      }
      printf("\n");
      
    }
    
  } else {
    
    len = strlen(out_file);
    name = ARRAY_ALLOC((len + 1), char);
    strcpy(name, out_file);
    name[len] = '\0';
    
    if ((out = fopen(name, "w")) == NULL) {
      fprintf(stderr, "Error in writing file %s\n", out_file);
      exit(1);
    }
    
    for (i = 0; i < res->a; i++) {
      fprintf(out, "\t%s", res->attributes[i]);
    }
    fprintf(out, "\n");
    
    for (i = 0; i < res->n; i++) {
      
      idx = perm[i];
      fprintf(out, "%s", res->variables[idx]);
      
      v = res->vectors[idx];
      for (j = 0; j < res->a; j++) {
	fprintf(out, "\t%e", v->data[j]);
      }
      fprintf(out, "\n");
      
    }
    
    fclose(out);
    free(name);
    
  }
  
}

void
print_PCA_results(char *out_file,
		  PCA_RESULT *res,
		  int axis)
{
  
  char *name;
  RESULT *r;

  if (out_file == NULL) {

    if (axis == 1 || axis == 2) {

      if (axis == 2) {
	printf("Samples\n");
      }

      r = res->sres;
      print_results(out_file, r, 0, 0);

    }

    if (axis == 0 || axis == 2) {  

      if (axis == 2) {
	printf("Variables\n");
      }

      r = res->gres;
      print_results(out_file, r, 0, 0);

    }

  } else {

    name = ARRAY_ALLOC((strlen(out_file) + strlen(".xxx") + 1), char);

    if (axis == 1 || axis == 2) {

      r = res->sres;
      strcpy(name, out_file);
      strcat(name, ".smp");
      print_results(name, r, 0, 0);

    }

    if (axis == 0 || axis == 2) {  

      r = res->gres;
      strcpy(name, out_file);
      strcat(name, ".gne");
      print_results(name, r, 0, 0);

    }

    free(name);

  }

}

void
print_mult_results(int n_res,
		   char **out_file,
		   RESULT **res,
		   int order,
		   int idx)
{

  int i;
  
  for (i=0; i<n_res; i++) {
    print_results(out_file[i],
		  res[i],
		  order,
		  idx);
  }

}

int 
count_lines(char *file_name)
{

  FILE *in;
  int i = 0;
  char line[8000000];
  
  if ((in = fopen(file_name, "r")) == NULL) { 
    fprintf(stderr, "Error in opening file %s\n", file_name);
    perror(0);
    exit(1);
  }
  
  while (fgets(line, sizeof(line), in)) {
    i++;
  }
  
  fclose(in);
  
  return i;
  
}

int 
count_wells_ks(char *file_name)
{

  FILE *in;
  int i, idx, tlen, w = 0;
  char line[8000000];
  char *tok;
  char **wells;

  wells = ARRAY_ALLOC(5000, char *);
  
  if ((in = fopen(file_name, "r")) == NULL) { 
    fprintf(stderr, "Error in opening file %s\n", file_name);
    perror(0);
    exit(1);
  }
  
  if (fgets(line, sizeof(line), in) == NULL) {
    fprintf(stderr, "Premature end of file %s\n", file_name);
    exit(1);
  }
  
  while (fgets(line, sizeof(line), in)) {
    idx = -1;
    tok = strtok (line, "\t");
    tok = strtok (tok, ".");
    if (w > 5000) {
      fprintf(stderr, "Too many wells in file %s\n", file_name);
      exit(1);
    } else {
      if (w > 0) {
	for (i = 0; i < w; i++) {
	  if (!strcmp(tok, wells[i])) {
	    idx = i;
	    break;
	  }
	}
	if (idx < 0) {
	  tlen = strlen(tok);
	  wells[w] = ARRAY_ALLOC((tlen+1), char);
	  strncpy(wells[w], tok, tlen);
	  wells[w][tlen] = '\0';
	  w++;
	}
      } else {
	tlen = strlen(tok);
	wells[w] = ARRAY_ALLOC((tlen+1), char);
	strncpy(wells[w], tok, tlen);
	wells[w][tlen] = '\0';
	w++;
      }
    }  
  }

  fclose(in);

  free(wells);

  return w;
  
}

int 
count_fields(char *file_name)
{

  FILE *in;
  int i = 0;
  char line[8000000];
  char *tok;
  
  if ((in = fopen(file_name, "r")) == NULL) { 
    fprintf(stderr, "Error in opening file %s\n", file_name);
    perror(0);
    exit(1);
  }
  
  if (fgets(line, sizeof(line), in) == NULL) {
    fprintf(stderr, "Premature end of file %s\n", file_name);
    exit(1);
  }
  
  tok = strtok (line, "\t");
  while (tok != NULL) {
    i++;
    tok = strtok (NULL, "\t\n");
  }
  
  fclose(in);
  
  return i;
  
}

int 
count_not_null_fields(char *file_name)
{

  FILE *in;
  int i = 0;
  char line[8000000];
  char *tok;
  
  // This should be a factor file
  if ((in = fopen(file_name, "r")) == NULL) { 
    fprintf(stderr, "Error in opening file %s\n", file_name);
    perror(0);
    exit(1);
  }
  
  // First line is the header
  if (fgets(line, sizeof(line), in) == NULL) {
    fprintf(stderr, "Premature end of file %s\n", file_name);
    exit(1);
  }
  
  // The second line is the first factor
  if (fgets(line, sizeof(line), in) == NULL) {
    fprintf(stderr, "Premature end of file %s\n", file_name);
    exit(1);
  }

  // The first token is the name of the factor,
  // the rest of the tokens are the data we
  // need to count
  tok = strtok (line, "\t");
  tok = strtok (NULL, "\t\n");
  while (tok != NULL) {
    if (strcmp(tok, "Unassigned") == 0) {
      i++;
    }
    tok = strtok (NULL, "\t\n");
  }
  
  fclose(in);
  
  return i;
  
}

KSDATA *
read_ksdata(char *file_name,
	    int n_features,
            int n_wells)
{
  FILE *in;
  int i, j, idx, c, w, tlen;
  int max = 8000000; // maximum number of cells per well
  double **tmp;
  KSDATA *ks;
  char line[8000000];
  char row[18];
  char *tok;

  ks = TYPE_ALLOC(KSDATA);
  ks->w = n_wells;
  ks->f = n_features;
  ks->c = ARRAY_ALLOC(ks->w, int);
  ks->wells = ARRAY_ALLOC(ks->w, char *);
  ks->features = ARRAY_ALLOC(ks->f, char *);
  ks->data = ARRAY_ALLOC(ks->f, double **);
  for (i = 0; i < ks->f; i++) {
    ks->data[i] = ARRAY_ALLOC(ks->w, double *);
  }

  if ((in = fopen(file_name, "r")) == NULL) { 
    fprintf(stderr, "Error in opening file %s\n", file_name);
    perror(0);
    exit(1);
  }
  
  /* Read the header */
  
  if (fgets(line, sizeof(line), in) == NULL) {
    fprintf(stderr, "Premature end of file %s\n", file_name);
    exit(1);
  }
  
  i = 0;
  tok = strtok (line, "\t\n");
  while (tok != NULL) {
    tlen = strlen(tok);
    ks->features[i] = ARRAY_ALLOC((tlen+1), char);
    strncpy(ks->features[i], tok, tlen);
    ks->features[i][tlen] = '\0';
    i++;
    tok = strtok (NULL, "\t\n");
  }
  if (i != n_features) {
    fprintf(stderr, "Bad feature name count in file %s\n", file_name);
    exit(1);
  }
  
  /* Read the rest of the file */
  
  // Temporary array to put the data for one well and many features
  tmp = ARRAY_ALLOC(ks->f, double *);
  for (i = 0; i < ks->f; i++) {
    tmp[i] = ARRAY_ALLOC(max, double);
  }

  // We start at well 0;
  w = 0;
  while (fgets(line, sizeof(line), in)) {
    idx = -1;
    // Get the well identifier for the row 
    tok  = strtok (line, ".");
    if (w > n_wells ) {
      fprintf(stderr, "There are more wells than specified in file %s\n", file_name);
      exit(1);
    } else {
      if (w > 0) {
	// A previously defined well
	for (i = 0; i < w; i++) {
	  if (!strcmp(tok, ks->wells[i])) {
	    idx = i;
	    break;
	  }
	}
	// A new well
	if (idx < 0) {
	  // Copy the data to the ks object
	  for (i = 0; i < ks->f; i++) {
	    ks->data[i][w-1] = ARRAY_ALLOC(c, double);
	    for (j = 0; j < c; j++) {
	      ks->data[i][w-1][j] = tmp[i][j];
	    }
	  }
	  // Number of cells in the well
	  ks->c[w-1] = c;
	  // Set the new well identifier
	  tlen = strlen(tok);
	  ks->wells[w] = ARRAY_ALLOC((tlen+1), char);
	  strncpy(ks->wells[w], tok, tlen);
	  ks->wells[w][tlen] = '\0';
	  idx = w;
	  c = 0;
	  w++;
	}
      } else {
	// The very first line of the file
	// We set the well identifier
	tlen = strlen(tok);
	ks->wells[w] = ARRAY_ALLOC((tlen+1), char);
	strncpy(ks->wells[w], tok, tlen);
	ks->wells[w][tlen] = '\0';
	idx = w;
	c = 0;
	w++;
      }
    }
    // The reset of the line
    strcpy(row, tok);
    strcat(row, ".");
    tok = strtok (NULL, "\t\n");
    strcat(row, tok);
    strcat(row, ".");
    tok = strtok (NULL, "\t\n");
    for (i = 0; i < n_features; i++) {
      if (tok != NULL) {
	tmp[i][c] = (double)atof(tok);
      } else {
	tmp[i][c] = MD;
      }
      tok = strtok (NULL, "\t\n");
    }
    c++;
  }
  
  // Copy the last well into the ks object
  for (i = 0; i < ks->f; i++) {
    ks->data[i][w-1] = ARRAY_ALLOC(c, double);
    for (j = 0; j < c; j++) {
      ks->data[i][w-1][j] = tmp[i][j];
    }
  }
  // Number of cells in the last well
  ks->c[w-1] = c;

  fclose(in);

  free(tmp);

  return ks;
    
}
