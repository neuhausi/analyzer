/*********************************************************************
 Copyright 2006 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: svd.c,v 1.4 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#define EXTERN
#include "svd.h"

void 
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **factor_file,
		char **output_file,
		int *samples,
		int *variables,
		int *components,
		int *eigenvs,
		enum METHOD_ENUM *method,
		int *transpose)
{

  int c;
  extern char *optarg;
  int errflg = 0;

  *data_file = NULL;

  *factor_file = NULL;

  *output_file = NULL;

  *samples = 0;

  *variables = 0;

  *components = 3;

  *eigenvs = 0;

  *method = COR;

  *transpose = 0;

  while ((c = getopt(argc, argv, "d:f:o:s:g:c:em:tvD")) != EOF)
    switch (c) {
    case 'd':
      *data_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*data_file, optarg);
      break;
    case 'f':
      *factor_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*factor_file, optarg);
      break;
    case 'o':
      *output_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*output_file, optarg);
      break;
    case 's':
      *samples = atoi(optarg);
      break;
    case 'g':
      *variables = atoi(optarg);
      break;
    case 'c':
      *components = atoi(optarg);
      break;
    case 'e':
      *eigenvs = 1;
      break;
    case 'm':
      if (strcmp(optarg, "cov") == 0) {
	*method = COV;
      } else {
	*method = COR;
      }
      break;
    case 't':
      *transpose = 1;
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

  if (*data_file == NULL) { 
    errflg++;
  }

  if (errflg) {
    printf("\nUsage: svd <options> (options in square brackets are optional)\n\n");
    printf(" -d  <string>   data file. A tab delimeted file with the following format.\n");
    printf("                A header line with sample names. Lines with expression values\n");
    printf("                preceded by the variable name at the begining of the line.\n");
    printf(" -f [<string>]  factor file. A tab delimet file containing the data to include.\n");
    printf("                If not specified, all the data is included.\n");
    printf(" -o [<string>]  output file. If no file is specified it defaults to stdout.\n");
    printf("                Two output files will be generated one with the extension '.gne'.\n");
    printf("                for the variables and one with extension '.smp' for the samples.\n");
    printf(" -s [<integer>] number of samples. It should correspond to the number of names\n");
    printf("                in the first line of the data file.\n");
    printf(" -g [<integer>] number of variables. It should correspond to the number of lines + 1\n");
    printf("                in the data file.\n");
    printf(" -c [<integer>] number of components (defaults to 3).\n");
    printf(" -e [<switch>]  output eigenvalue, variance explained and eigenvectors.\n");
    printf(" -m [<string>]  [cor|cov]. Type of pre-processing. Correlation or Covariance\n");
    printf("                (defaults to cor).\n");
    printf(" -t [<switch>]  transpose data matrix.\n\n");
    printf("Examples\n");
    printf("\n\tsvd -d ./t/svd.dat\n\n");
    exit(0);
  }
   
}

RESULT *
run_svd(MATRIX *matrix,
	int components,
	int eigenvs,
	enum METHOD_ENUM method)
{

  int i, j, k, n, m, max, ac, len;
  double **v, *w, sum;
  RESULT *res;
  VECTOR *vd;	

  m = matrix->m;
  n = matrix->n;

  // Memory Allocation

  res = set_res_object(m, n, components, eigenvs);

  v  = ARRAY_ALLOC(m, double *);
  w  = ARRAY_ALLOC(m, double);

  for (i=0;i<m;i++) {
    v[i] = ARRAY_ALLOC(m, double);
  }

  if (method == COR) {
    corcol(matrix->data, n, m);
  } else if (method == COV) {
    covcol(matrix->data, n, m);
  }

  svdcmp(matrix->data,n,m,w,v);

  max = components > m ? m : components;

  if (eigenvs == 1) {

    // add the Eigenvalues
    sum=0.0;
    for (i=0;i<m;i++)  {
      sum += w[i];
    }

    // Eigenvalues and Variance Explained
    k=0;
    for (i=0;i<m;i++) {
      vd = res->vectors[i];
      vd->data[0] = w[i];
      vd->data[1] = (w[i] / sum) * 100;      
    }

    // Eigenvectors
    for (i=0;i<m;i++) {
      ac=2;
      vd = res->vectors[i];
      for (j=0;j<max;j++)  {
	vd->data[ac] = v[i][j] * w[j];
	ac++;
      }
    }

  }

  for (i=0;i<m;i++) {
    // Set the sample names
    len = strlen(matrix->cols[i]);
    res->variables[i] = ARRAY_ALLOC((len+1), char);
    strncpy(res->variables[i], matrix->cols[i], len);
    res->variables[i][len] = '\0';
    // Set the data
    vd = res->vectors[i];
    ac = eigenvs == 1 ? max+2 : 0;
    for (j=0;j<max;j++)  {
      vd->data[ac] = v[i][j];
      ac++;
    }
  }

  free(v);
  free(w);

  return res;
 
}

RESULT *
set_res_object(int s,
	       int n,
	       int c,
	       int e)
{

  int max, ea, len, cnt, i;
  RESULT *res;
  VECTOR *v;
  char strc[8];
  char str[18];

  res = TYPE_ALLOC(RESULT);

  max = c > s ? s : c;
  ea = e == 1 ? (max*2) + 2 : max;
  res->n = s;
  res->a = ea;
  res->attributes = ARRAY_ALLOC(res->a, char *);
  res->variables = ARRAY_ALLOC(s, char *);
  res->vectors = ARRAY_ALLOC(s, VECTOR *);
  for (i=0;i<s;i++) {
    v = res->vectors[i] = TYPE_ALLOC(VECTOR);
    v->data = ARRAY_ALLOC(ea, double);
    v->n = ea;
  }
  cnt = 0;
  if (e == 1) {
    strcpy(str, "Eigenvalue");
    len = strlen(str);
    res->attributes[cnt] = ARRAY_ALLOC((len+1), char);
    strncpy(res->attributes[cnt], str, len);
    res->attributes[cnt][len] = '\0';
    cnt++;
    strcpy(str, "%VarianceExplained");
    len = strlen(str);
    res->attributes[cnt] = ARRAY_ALLOC((len+1), char);
    strncpy(res->attributes[cnt], str, len);
    res->attributes[cnt][len] = '\0';
    cnt++;
    for (i=0;i<max;i++) {
      strcpy(str, "Eigenvector");
      sprintf(strc, "%i", i+1);
      strcat(str, strc);
      len = strlen(str);
      res->attributes[cnt] = ARRAY_ALLOC((len+1), char);
      strncpy(res->attributes[cnt], str, len);
      res->attributes[cnt][len] = '\0';
      cnt++;
    }
  }
  for (i=0;i<max;i++) {
    strcpy(str, "PC");
    sprintf(strc, "%i", i+1);
    strcat(str, strc);
    len = strlen(str);
    res->attributes[cnt] = ARRAY_ALLOC((len+1), char);
    strncpy(res->attributes[cnt], str, len);
    res->attributes[cnt][len] = '\0';
    cnt++;
  }
  
  return res;

}

void
corcol (double **data, int n, int m)
{
  
  // This code was adapted from a public domain version by
  // F. Murtagh and A. Heck, Multivariate Data Analysis, Kluwer
  // Academic, Dordrecht, 1987 

  double x, *mean, *stddev;
  int i, j;
  
  mean = ARRAY_ALLOC(m, double);
  stddev = ARRAY_ALLOC(m, double);

  // Determine mean of column vectors of input data matrix

  for (j=0;j<m;j++) {
    mean[j]=0.0;
    for (i=0;i<n;i++) {
      mean[j] += data[i][j];
    }
    mean[j] /= (double)n;
  }

  // Determine standard deviations of column vectors of data matrix

  for (j=0;j<m;j++){
    stddev[j] = 0.0;
    for (i=0;i<n;i++){
      stddev[j] += ((data[i][j]-mean[j])*(data[i][j]-mean[j]));
    }
    stddev[j] /= (double)n;
    stddev[j] = sqrt(stddev[j]);
    if (stddev[j] <= TINY) stddev[j]=1.0;
  }

  // Center and reduce the column vectors

  for (i=0;i<n; i++){
    for (j=0;j<m;j++){
      data[i][j] -= mean[j];
      x = sqrt((double)n);
      x *= stddev[j];
      data[i][j] /= x;
    }
  }

}

void
covcol (double **data, int n, int m)
{

  // This code was adapted from a public domain version by
  // F. Murtagh and A. Heck, Multivariate Data Analysis, Kluwer
  // Academic, Dordrecht, 1987 

  double *mean;
  int i, j;
   
  mean = ARRAY_ALLOC(m, double);

  // Determine mean of column vectors of input data matrix

  for (j=0;j<m;j++) {
    mean[j]=0.0;
    for (i=0;i<n;i++) {
      mean[j] += data[i][j];
    }
    mean[j] /= (double)n;
  }
  
  // Center the column vectors
  
  for (i=0;i<n;i++){
    for (j=0;j<m;j++){
      data[i][j] -= mean[j];
    }
  }
  
}

int
main(int argc,
     char **argv)

{

    char *data_file, *factor_file, *output_file;
    int  samples, variables, components, eigenvs, transpose, binary;
    enum METHOD_ENUM method;
    MATRIX *matrix;
    FACTOR *factor;
    RESULT *res;

    verbose = 0;

    debug = 0;

    binary = 0;

    parse_arguments(argc, argv, &data_file, &factor_file, &output_file, 
		    &samples, &variables, &components, &eigenvs,
		    &method, &transpose);

    binary = is_binary_file(data_file);

    if (binary > 0) {
      matrix = read_binary_matrix(data_file);
    } else {
      if (samples == 0) {
 	samples = count_fields(data_file);
      }
      if (variables == 0) {
 	variables = count_lines(data_file);
	variables--;
      }
      matrix = read_matrix(data_file, samples, variables);
    }

    if (transpose == 1) {
      matrix = transpose_matrix(matrix);
    }

    if (factor_file != NULL) {
      factor = read_factor(factor_file, 1, matrix->m);
      matrix = remove_null_matrix(matrix, factor);
    }

    res = run_svd(matrix, components, eigenvs, method);

    print_results(output_file, res, 0, 0); 
    
    free(data_file);
    free(factor_file);
    free(output_file);
    free(matrix);
    free(res);

    exit(0);

}
