/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: correct.c,v 1.8 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#define EXTERN
#include "correct.h"

void
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **output_file,
		int *cols,
		int *variables,
		int **inds,
		double *Q,
		enum TT *tt,
		int *n_inds)
{
  
  int c, max_inds;
  extern char *optarg;
  int errflg = 0;
  
  *data_file = NULL;
  
  *output_file = NULL;
  
  /* Maximum of 10 columns to calculate */
  max_inds = 10;

  *inds = ARRAY_ALLOC(max_inds, int);

  *cols = 0;

  *variables = 0;

  *Q = 0.05;
  
  *tt = FDR;

  *n_inds = 0;
  
  while ((c = getopt(argc, argv, "d:o:i:g:c:t:q:h:vD")) != EOF)
    switch (c) {
    case 'd':
      *data_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*data_file, optarg);
      break;
    case 'o':
      *output_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*output_file, optarg);
      break;
    case 'i':
      if (*n_inds <= max_inds) {
	(*inds)[*n_inds] = atoi(optarg);
	*n_inds +=1;
      }
      break;
    case 'g':
      *variables = atoi(optarg);
      break;
    case 'c':
      *cols = atoi(optarg);
      break;
    case 't':
      if (strcmp(optarg, "bonferroni") == 0) {
	*tt = BONFERRONI;
      } else {
	*tt = FDR;
      }
      break;
    case 'q':
      *Q = atof(optarg);
      break;
    case 'h':
      errflg++;
      break;
    case 'v':
      verbose=1;
      break;
    case 'D':
      debug=1;
      break;
    case '?':
      errflg++;
      break;
    }

  if (*data_file == NULL || *n_inds == 0) { 
    errflg++;
  }

  if (errflg) {
    printf("\nUsage: correct <options>  (options in square brackets are optional)\n\n");
    printf(" -d  <string>    data file. A tab delimeted file with the following format.\n");
    printf("                 A header line with column names. Lines with expression values\n");
    printf("                 preceded by the variable name at the begining of the line.\n");
    printf(" -o  [<string>]  output file. If no file is specified it defaults to stdout.\n");
    printf(" -g  [<integer>] number of variables. It should correspond to the number of lines + 1\n");
    printf("                 in the data file.\n");
    printf(" -c  [<integer>] number of columns. It should correspond to the number of names\n");
    printf("                 in the first line of the data file.\n");
    printf(" -i  <integer>   indeces of columns to use for the calculations (columns start at 0)\n");
    printf(" -t  [<string>]  [fdr|bonferroni]. Type of correction. (default is fdr).\n");
    printf(" -q  [<float>]   q value to use in false discovery rate. (default is 0.05).\n\n");
    printf("Examples\n");
    printf("\n\tcorrect -d ./t/correct.dat -i 2\n\n");
    exit(0);
  }
   
}

RESULT *
correct(MATRIX *matrix,
	double Q,
	enum TT tt, 
	int *inds,
	int n_cols,
	double **cutoff,
	int  **nsig)
{

  int i, j, len;
  double cdata[n_cols][matrix->n];
  char str[1000], strc[10];
  RESULT *res;
  VECTOR *vr;	

  res = TYPE_ALLOC(RESULT);
  res->n = matrix->n;
  res->a = matrix->m + n_cols;
  res->attributes = ARRAY_ALLOC(res->a, char *);
  res->variables = ARRAY_ALLOC(res->n, char *);
  res->vectors = ARRAY_ALLOC(res->n, VECTOR *);

  *cutoff = ARRAY_ALLOC(n_cols, double);
  *nsig = ARRAY_ALLOC(n_cols, int);

  for (i=0;i<n_cols;i++) {
    strcpy(str, matrix->cols[inds[i]]);
    if (tt == FDR) {
      strcat(str, " FDR (");
      sprintf(strc, "%4.2f", Q);
      strcat(str, strc);
      strcat(str, ")");
    } else {
      strcat(str, " bonferroni");
    }
    len = strlen(str);
    res->attributes[i] = ARRAY_ALLOC((len+1), char);
    strncpy(res->attributes[i], str, len);
    res->attributes[i][len] = '\0';
  }
  for (i=0;i<matrix->m;i++) {
    strcpy(str, matrix->cols[i]);
    len = strlen(str);
    res->attributes[i+n_cols] = ARRAY_ALLOC((len+1), char);
    strncpy(res->attributes[i+n_cols], str, len);
    res->attributes[i+n_cols][len] = '\0';
  }   
  
  matrix = transpose_matrix(matrix);
  
  for (i=0;i<n_cols;i++) {
    if (tt == FDR) {
      fdr(res->n, matrix->data[inds[i]], cdata[i], Q, 1, &(*cutoff)[i], &(*nsig)[i]);
      if (verbose || debug) {
	strcpy(str, matrix->rows[inds[i]]);
	strcat(str, " FDR (");
	sprintf(strc, "%5.3f", Q);
	strcat(str, strc);
	strcat(str, ")");
	printf("%s\n# p-values     Cutoff   # signif\n", str);
	printf("%10i %10.3g %10i\n", res->n, (*cutoff)[i], (*nsig)[i]);
      }
    } else {
      bonferroni(matrix->data[inds[i]], cdata[i], res->n);
    }
  }

  matrix = transpose_matrix(matrix);

  for (i = 0; i < res->n; i++) {

    len = strlen(matrix->rows[i]);
    res->variables[i] = ARRAY_ALLOC((len+1), char);
    strncpy(res->variables[i], matrix->rows[i], len);
    res->variables[i][len] = '\0';

    vr = res->vectors[i] = TYPE_ALLOC(VECTOR);
    vr->data = ARRAY_ALLOC(res->a, double);
    vr->n = res->a;

    for (j=0;j<n_cols;j++) {
      vr->data[j] = cdata[j][i];
    }

    for (j=0;j<matrix->m;j++) {
      vr->data[j+n_cols] = matrix->data[i][j];
    }

  }

  return res;

}

void
print_FDR_summary(char *out_file,
		  RESULT *res,
		  double *cutoff,
		  int *nsig,
		  int n_cols)
{

  FILE *out;
  char *name;
  int i;

  if (out_file == NULL) {

    for (i=0;i<n_cols;i++) {
      printf("\t%s", res->attributes[i]);
    }
    printf("\n");

    printf("p-values");
    for (i=0;i<n_cols;i++) {
      printf("\t%i", res->n);
    }
    printf("\nCutoff");
    for (i=0;i<n_cols;i++) {
      printf("\t%e", cutoff[i]);
    }

    printf("\n# signif");
    for (i=0;i<n_cols;i++) {
      printf("\t%i", nsig[i]);
    }
    printf("\n");


  } else {

    name = ARRAY_ALLOC((strlen(out_file) + strlen(".sum") + 1), char);
    strcpy(name, out_file);
    strcat(name, ".sum");
    
    if ((out = fopen(name, "w")) == NULL) {
      fprintf(stderr, "Error in writing file %s\n", out_file);
      perror(0);
      exit(1);
    }
 
    for (i=0;i<n_cols;i++) {
      fprintf(out, "\t%s", res->attributes[i]);
    }
    fprintf(out, "\n");

    fprintf(out, "p-values");
    for (i=0;i<n_cols;i++) {
      fprintf(out, "\t%i", res->n);
    }
    fprintf(out, "\nCutoff");
    for (i=0;i<n_cols;i++) {
      fprintf(out, "\t%e", cutoff[i]);
    }

    fprintf(out, "\n# signif");
    for (i=0;i<n_cols;i++) {
      fprintf(out, "\t%i", nsig[i]);
    }
    fprintf(out, "\n");

    fclose(out);
    free(name);

  }
   


}

int
main(int argc,
     char **argv)

{

  char *data_file, *output_file;
  int  cols, variables, n_cols, *inds, *nsig;
  double Q, *cutoff;
  enum TT tt;
  MATRIX *matrix;
  RESULT *res;
  
  verbose = 0;

  debug = 0;

  parse_arguments(argc, argv, &data_file, &output_file,
		  &cols, &variables, &inds, &Q, &tt, &n_cols);
  
  if (cols == 0) {
    cols = count_fields(data_file);
  }
  
  if (variables == 0) {
    variables = count_lines(data_file);
    variables--;
  }
  
  matrix = read_matrix(data_file, cols, variables);

  res = correct(matrix, Q, tt, inds, n_cols, &cutoff, &nsig);

  print_results(output_file, res, 1, inds[0]+n_cols);

  if (tt == FDR) {
    print_FDR_summary(output_file, res, cutoff, nsig, n_cols);
  }

  exit(0);

}
