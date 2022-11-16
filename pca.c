/*********************************************************************
 Copyright 2004 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Principal Components Analysis or the Karhunen-Loeve expansion is a
   classical method for dimensionality reduction or exploratory data
   analysis.  This code was adapted from a public domain version by
   F. Murtagh and A. Heck, Multivariate Data Analysis, Kluwer
   Academic, Dordrecht, 1987 . Also several subrutines were also
   adapted from Numerical Recipes in C, W H Press et al., 1988

 $Id: pca.c,v 1.17 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#define EXTERN
#include "pca.h"

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
		enum AXIS *axis,
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

  *axis = BOTH;

  *transpose = 0;

  while ((c = getopt(argc, argv, "d:f:o:s:g:c:em:a:tvD")) != EOF)
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
      } else if (strcmp(optarg, "sscp") == 0) {
	*method = SSCP;
      } else {
	*method = COR;
      }
      break;
    case 'a':
      if (strcmp(optarg, "variables") == 0) {
	*axis = VARIABLES;
      } else if (strcmp(optarg, "samples") == 0) {
	*axis = SAMPLES;
      } else {
	*axis = BOTH;
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
    printf("\nUsage: pca <options> (options in square brackets are optional)\n\n");
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
    printf(" -m [<string>]  [cor|cov|sscp]. Type of PCA. Correlation, Covariance or Sum of\n");
    printf("                Squares / Cross Products (defaults to cor).\n");
    printf(" -a [<string>]  [variables|samples|both]. Axis to run the PCA. Variables, samples or both\n");
    printf("                (defaults to both).\n");
    printf(" -t [<switch>]  transpose data matrix.\n\n");
    printf("Examples\n");
    printf("\n\tpca -d ./t/pca.dat\n\n");
    exit(0);
  }
   
}

PCA_RESULT *
run_pca(MATRIX *matrix,
	int components,
	int eigenvs,
	enum METHOD_ENUM method,
	enum AXIS axis)
{

  int i, j, k, l, m, n, max, ac, len, *idx;
  double **symmat, **symmat2, *evals, *interm, sum;
  PCA_RESULT *res, *set_res_object(int s, int n, int c, int e);
  RESULT *sr;
  RESULT *gr;
  VECTOR *v;	

  m = matrix->m;
  n = matrix->n;

  // Memory Allocation

  res = set_res_object(m, n, components, eigenvs);
  sr = res->sres;
  gr = res->gres;

  if (axis == VARIABLES || axis == BOTH) {  
    for (i=0;i<n;i++) {
      // Set the variable names
      len = strlen(matrix->rows[i]);
      gr->variables[i] = ARRAY_ALLOC((len+1), char);
      strncpy(gr->variables[i], matrix->rows[i], len);
      gr->variables[i][len] = '\0';
    }
  }

  if (axis == SAMPLES || axis == BOTH) {
    for (j=0;j<m;j++) {
      // Set the sample names
      len = strlen(matrix->cols[j]);
      sr->variables[j] = ARRAY_ALLOC((len+1), char);
      strncpy(sr->variables[j], matrix->cols[j], len);
      sr->variables[j][len] = '\0';
    }
  }

  idx = ARRAY_ALLOC(m, int);
  symmat = ARRAY_ALLOC(m, double *);
  symmat2 = ARRAY_ALLOC(m, double *);
  evals = ARRAY_ALLOC(m, double);
  interm = ARRAY_ALLOC(m, double);
  for (i=0;i<m;i++) {
    symmat[i] = ARRAY_ALLOC(m, double);
    symmat2[i] = ARRAY_ALLOC(m, double); 
  }

  if (method == COV) {
    covcol(matrix->data, n, m, symmat);
  } else if (method == SSCP) {
    scpcol(matrix->data, n, m, symmat);
  } else {
    corcol(matrix->data, n, m, symmat);
  }

  if (axis == SAMPLES || axis == BOTH) {
    for (i=0;i<m;i++) {
      for (j=0;j<m;j++) {
	symmat2[i][j]=symmat[i][j];
      }
    }
  }

  // Triangular decomposition
  tred2(symmat, m, evals, interm);

  // Reduction of sym. trid. matrix
  tqli(evals, interm, m, symmat);
 
  sort_by_index(m, evals, idx);

  max = components > m ? m : components;

  if (eigenvs == 1) {

    // add the Eigenvalues
    sum=0.0;
    for (i=0;i<m;i++)  {
      sum += evals[i];
    }

    // Eigenvalues and Variance Explained
    k=0;
    for (j=(m-1);j>=0;j--) {
      v = sr->vectors[k];
      v->data[0] = evals[idx[j]];
      v->data[1] = (evals[idx[j]] / sum) * 100;      
      k++;
    }

    // Eigenvectors
    for (j=0;j<m;j++) {
      ac=2;
      v = sr->vectors[j];
      for (i=1;i<=max;i++)  {
	v->data[ac] = symmat[j][idx[m-i]];
	ac++;
      }
    }

  }

  if (axis == VARIABLES || axis == BOTH) {  

    max = components > n ? n : components;

    // Form projections of variables on principal components
    // Store in 'matrix->data', overwriting original data.
    for (i=0;i<n;i++) {
      for (j=0;j<m;j++) {
	if (! isnan(matrix->data[i][j])) {	
	  interm[j] = matrix->data[i][j];
	}
      }
      for (k=0;k<=m;k++) {
	matrix->data[i][k]=0.0;
	for (l=0;l<m;l++) {
	  if (! isnan(matrix->data[i][k])) {	
	    matrix->data[i][k] += interm[l] * symmat[l][m-k];
	  }
	}
      }
    }

    // Projections of variables on principal components
    for (i=0;i<n;i++) {
      // Set the data
      v = gr->vectors[i];
      ac=0;
      for (j=1;j<=max;j++)  {
	if (! isnan(matrix->data[i][m-idx[m-j]])) {	
	  v->data[ac] = matrix->data[i][m-idx[m-j]];
	  ac++;
	}
      }
    }
  }

  if (axis == SAMPLES || axis == BOTH) {

    max = components > m ? m : components;

    // Form projections of samples on principal components
    // Store in 'symmat2', overwriting what was stored in this.
    for (j=0;j<m;j++) {
      for (k=0;k<m;k++) {
	interm[k]=symmat2[j][k];
      }
      for (i=0;i<=m;i++) {
	symmat2[j][i]=0.0;
	for (l=0;l<m;l++) {
	  symmat2[j][i] += interm[l] * symmat[l][m-i];
	}
	if (evals[m-i] > TINY) {
	  symmat2[j][i] /= sqrt(evals[m-i]);
	} else {
	  symmat2[j][i]=0.0;
	}
      }
    }
 
    // Projections of samples on principal components
    for (j=0;j<m;j++) {
      // Set the data
      v = sr->vectors[j];
      ac = eigenvs == 1 ? max+2 : 0;
      for (k=1;k<=max;k++)  {
	v->data[ac] = symmat2[j][m-idx[m-k]];
	ac++;
      }
    }

  }

  free(symmat);
  free(symmat2);
  free(evals);
  free(interm);
  free(idx);

  return res;
 
}

PCA_RESULT *
set_res_object(int s,
	       int n,
	       int c,
	       int e)
{

  int max, ea, len, cnt, i;
  PCA_RESULT *res;
  RESULT *sr;
  RESULT *gr;
  VECTOR *v;
  char strc[8];
  char str[18];

  res = TYPE_ALLOC(PCA_RESULT);

  sr = res->sres = TYPE_ALLOC(RESULT);
  gr = res->gres = TYPE_ALLOC(RESULT);

  // Setup the samples first

  max = c > s ? s : c;
  ea = e == 1 ? (max*2) + 2 : max;
  sr->n = s;
  sr->a = ea;
  sr->attributes = ARRAY_ALLOC(sr->a, char *);
  sr->variables = ARRAY_ALLOC(s, char *);
  sr->vectors = ARRAY_ALLOC(s, VECTOR *);
  for (i=0;i<s;i++) {
    v = sr->vectors[i] = TYPE_ALLOC(VECTOR);
    v->data = ARRAY_ALLOC(ea, double);
    v->n = ea;
  }
  cnt = 0;
  if (e == 1) {
    strcpy(str, "Eigenvalue");
    len = strlen(str);
    sr->attributes[cnt] = ARRAY_ALLOC((len+1), char);
    strncpy(sr->attributes[cnt], str, len);
    sr->attributes[cnt][len] = '\0';
    cnt++;
    strcpy(str, "%VarianceExplained");
    len = strlen(str);
    sr->attributes[cnt] = ARRAY_ALLOC((len+1), char);
    strncpy(sr->attributes[cnt], str, len);
    sr->attributes[cnt][len] = '\0';
    cnt++;
    for (i=0;i<max;i++) {
      strcpy(str, "Eigenvector");
      sprintf(strc, "%i", i+1);
      strcat(str, strc);
      len = strlen(str);
      sr->attributes[cnt] = ARRAY_ALLOC((len+1), char);
      strncpy(sr->attributes[cnt], str, len);
      sr->attributes[cnt][len] = '\0';
      cnt++;
    }
  }
  for (i=0;i<max;i++) {
    strcpy(str, "PC");
    sprintf(strc, "%i", i+1);
    strcat(str, strc);
    len = strlen(str);
    sr->attributes[cnt] = ARRAY_ALLOC((len+1), char);
    strncpy(sr->attributes[cnt], str, len);
    sr->attributes[cnt][len] = '\0';
    cnt++;
  }
  
  // Setup the variables

  max = c > n ? n : c;
  cnt = 0;
  gr->n = n;
  gr->a = max;
  gr->attributes = ARRAY_ALLOC(gr->a, char *);
  gr->variables = ARRAY_ALLOC(n, char *);
  gr->vectors = ARRAY_ALLOC(n, VECTOR *);
  for (i=0;i<n;i++) {
    v = gr->vectors[i] = TYPE_ALLOC(VECTOR);
    v->data = ARRAY_ALLOC(max, double);
    v->n = max;
  }
  for (i=0;i<max;i++) {
    strcpy(str, "PC");
    sprintf(strc, "%i", i+1);
    strcat(str, strc);
    len = strlen(str);
    gr->attributes[cnt] = ARRAY_ALLOC((len+1), char);
    strncpy(gr->attributes[cnt], str, len);
    gr->attributes[cnt][len] = '\0';
    cnt++;
  }
  
  return res;

}

void
corcol (double **data, int n, int m, double **symmat)
{
  
  // This code was adapted from a public domain version by
  // F. Murtagh and A. Heck, Multivariate Data Analysis, Kluwer
  // Academic, Dordrecht, 1987 

  double x, mean[m], stddev[m];
  int i, j, j1, j2, ns[m];
  
  // Determine mean of column vectors of input data matrix

  for (j=0;j<m;j++) {
    mean[j]=0.0;
    ns[j]=0;
    for (i=0;i<n;i++) {
      if (! isnan(data[i][j])) {	
	mean[j] += data[i][j];
	ns[j]++;
      }
    }
    mean[j] /= (double)ns[j];
  }

  // Determine standard deviations of column vectors of data matrix

  for (j=0;j<m;j++){
    stddev[j] = 0.0;
    for (i=0;i<n;i++){
      if (! isnan(data[i][j])) {	
	stddev[j] += ((data[i][j]-mean[j])*(data[i][j]-mean[j]));
      }
    }
    stddev[j] /= (double)ns[j];
    stddev[j] = sqrt(stddev[j]);
    if (stddev[j] <= TINY) stddev[j]=1.0;
  }

  // Center and reduce the column vectors

  for (i=0;i<n; i++){
    for (j=0;j<m;j++){
      if (! isnan(data[i][j])) {	
	data[i][j] -= mean[j];
	x = sqrt((double)n);
	x *= stddev[j];
	data[i][j] /= x;
      }
    }
  }

  // Calculate the m * m correlation matrix

  for (j1=0;j1<m;j1++){
    symmat[j1][j1]=1.0;
    for (j2=j1+1;j2<m;j2++){
      symmat[j1][j2]=0.0;
      for (i=0;i<n;i++){
	if (! isnan(data[i][j1]) && ! isnan(data[i][j2])) {	
	  symmat[j1][j2] += (data[i][j1]*data[i][j2]);
	}
      }
      symmat[j2][j1] = symmat[j1][j2];
    }
  }

  if (debug) {
    printf ("\nCorrelation in columns\n");
    print_mat_d(symmat, m, m);
    printf ("\n");    
  }

}

void
covcol (double **data, int n, int m, double **symmat)
{

  // This code was adapted from a public domain version by
  // F. Murtagh and A. Heck, Multivariate Data Analysis, Kluwer
  // Academic, Dordrecht, 1987 

  double mean[m];
  int i, j, j1, j2, ns[m];
   
  // Determine mean of column vectors of input data matrix

  for (j=0;j<m;j++) {
    mean[j]=0.0;
    ns[j]=0;
    for (i=0;i<n;i++) {
      if (! isnan(data[i][j])) {	
	mean[j] += data[i][j];
	ns[j]++;
      }
    }
    mean[j] /= (double)ns[j];
  }
  
  // Center the column vectors
  
  for (i=0;i<n;i++){
    for (j=0;j<m;j++){
      if (! isnan(data[i][j])) {	
	data[i][j] -= mean[j];
      }
    }
  }
  
  // Calculate the m * m covariance matrix
  
  for (j1=0;j1<m;j1++){
    for (j2=j1;j2<m;j2++){
      symmat[j1][j2]=0.0;
      for (i=0;i<n;i++){
	if (! isnan(data[i][j1]) && ! isnan(data[i][j2])) {	
	  symmat[j1][j2] += data[i][j1]*data[i][j2];
	}
      }
      symmat[j2][j1]=symmat[j1][j2];
    }
  }
  
  if (debug) {
    printf ("\nCovariance in columns\n");
    print_mat_d(symmat, m, m);
    printf ("\n");    
  }

}

void
scpcol (double **data, int n, int m, double **symmat)
{

  // This code was adapted from a public domain version by
  // F. Murtagh and A. Heck, Multivariate Data Analysis, Kluwer
  // Academic, Dordrecht, 1987 

  int i, j1, j2;

  // Calculate the m * m sums-of-squares-and-cross-products matrix

  for (j1=0;j1<m;j1++){
    for (j2=j1;j2<m;j2++){
      symmat[j1][j2]=0.0;
      for (i=0;i<n;i++){
	if (! isnan(data[i][j1]) && ! isnan(data[i][j2])) {	
	  symmat[j1][j2] += data[i][j1]*data[i][j2];
	}
      }
      symmat[j2][j1]=symmat[j1][j2];
    }
  }
  
  if (debug) {
    printf ("\nDot product in columns\n");
    print_mat_d(symmat, m, m);
    printf ("\n");    
  }

}

int
main(int argc,
     char **argv)

{

    char *data_file, *factor_file, *output_file;
    int  samples, variables, components, eigenvs, transpose, binary, inta;
    enum METHOD_ENUM method;
    enum AXIS axis;
    MATRIX *matrix;
    FACTOR *factor;
    PCA_RESULT *res;

    verbose = 0;

    debug = 0;

    binary = 0;

    parse_arguments(argc, argv, &data_file, &factor_file, &output_file, 
		    &samples, &variables, &components, &eigenvs,
		    &method, &axis, &transpose);

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

    res = run_pca(matrix, components, eigenvs, method, axis);

    if (axis == VARIABLES) {
      inta = 0;
    } else if (axis == SAMPLES) {
      inta = 1;
    } else {
      inta = 2;
    }

    print_PCA_results(output_file, res, inta); 
    
    free(data_file);
    free(factor_file);
    free(output_file);
    free(matrix);
    free(res);

    exit(0);

}
