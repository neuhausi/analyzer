/*********************************************************************
 Copyright 2006 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

 $Id: fit.c,v 1.21 2009/09/17 21:30:01 neuhausi Exp $
**********************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#define EXTERN
#include "fit.h"

void
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **factor_file,
		char **covariate_file,
		char **output_file,
		int *samples,
		int *variables,
		char ***variable_names,
		int *variable_number,
		char **fit_factor,
		char **agg_factor,
		char **zero,
		int *logtrans,
		enum FUNCS *func,
		double **params,
		double *md,
		int *progress)
{
  
  int i, c, option_index, max_vars;
  extern char *optarg;
  int errflg = 0;
  
  *data_file = NULL;
  
  *factor_file = NULL;
  
  *covariate_file = NULL;

  *output_file = NULL;

  *samples = 0;
  
  *variables = 0;

  /* Maximum of 100 objects */
  max_vars = 100;

  *variable_names = ARRAY_ALLOC(max_vars, char *);
  
  *variable_number = 0;

  *fit_factor = NULL;
  
  *agg_factor = NULL;
  
  *zero = NULL;

  *logtrans = 0;

  *func = NLLS4P;

  *params = ARRAY_ALLOC(8, double);

  for (i = 0; i < 8; i++) {
    (*params)[i] = -999999;
  }

  *md = MD;
  
  *progress = 0;
  
  while (1) {

    static struct option long_options[] = {
      {"d", required_argument, 0, 'd'},
      {"f", required_argument, 0, 'f'},
      {"c", required_argument, 0, 'c'},
      {"o", required_argument, 0, 'o'},
      {"s", required_argument, 0, 's'},
      {"g", required_argument, 0, 'g'},
      {"n", required_argument, 0, 'n'},
      {"m", required_argument, 0, 'm'},
      {"a", required_argument, 0, 'a'},
      {"z", required_argument, 0, 'a'},
      {"t", required_argument, 0, 't'},
      {"e", required_argument, 0, 'e'},
      {"P1", required_argument, 0, 'P'},
      {"P2", required_argument, 0, 'Q'},
      {"P3", required_argument, 0, 'R'},
      {"P4", required_argument, 0, 'S'},
      {"P5", required_argument, 0, 'T'},
      {"P6", required_argument, 0, 'U'},
      {"P7", required_argument, 0, 'V'},
      {"P8", required_argument, 0, 'W'},
      {"p", no_argument, 0, 'p'},
      {"l", no_argument, 0, 'l'},
      {"v", no_argument, 0, 'v'},
      {"D", no_argument, 0, 'D'},
      {"h", no_argument, 0, 'h'},
      {0, 0, 0, 0}
    };

    option_index = 0;

    c = getopt_long_only(argc, argv, "d:f:c:o:s:g:n:m:a:z:t:P:Q:R:S:T:U:V:W:e:plvDh", long_options, &option_index);

    if (c == -1) break;

    switch (c) {

    case 'd':
      *data_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*data_file, optarg);
      break;
    case 'f':
      *factor_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*factor_file, optarg);
      break;
    case 'c':
      *covariate_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*covariate_file, optarg);
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
    case 'n':
      if (*variable_number <= max_vars) {
	(*variable_names)[*variable_number] = ARRAY_ALLOC((strlen(optarg) + 1), char);
	strcpy((*variable_names)[*variable_number], optarg);
	*variable_number += 1;
      }
      break;
    case 'm':
      *fit_factor = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*fit_factor, optarg);
      break;
    case 'a':
      *agg_factor = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*agg_factor, optarg);
      break;
    case 'z':
      *zero = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*zero, optarg);
      break;
    case 't':
      if (strcmp(optarg, "nlls4p") == 0) {
	*func = NLLS4P;
      } else {
	errflg++;
      }
      break;
    case 'l':
      *logtrans=1;
      break;
    case 'e':
      *md = atof(optarg);
      break;
    case 'P':
      (*params)[0] = atof(optarg);
      break;
    case 'Q':
      (*params)[1] = atof(optarg);
      break;
    case 'R':
      (*params)[2] = atof(optarg);
      break;
    case 'S':
      (*params)[3] = atof(optarg);
      break;
    case 'T':
      (*params)[4] = atof(optarg);
      break;
    case 'U':
      (*params)[5] = atof(optarg);
      break;
    case 'V':
      (*params)[6] = atof(optarg);
      break;
    case 'W':
      (*params)[7] = atof(optarg);
      break;
    case 'p':
      *progress=1;
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

  }

  if (*data_file == NULL || *factor_file == NULL || *covariate_file == NULL || 
      *fit_factor == NULL || *agg_factor == NULL) { 
    errflg++;
  }

  if (errflg) {
    printf("\nCurve fitting algorithm\n\n");
    printf("\nUsage: fit <options>  (options in square brackets are optional)\n\n");
    printf(" -d  <string>    data file. A tab delimeted file with the following format.\n");
    printf("                 A header line with sample names. Lines with expression values\n");
    printf("                 preceded by the variable name at the begining of the line.\n");
    printf(" -f  [<string>]  factor file. A tab delimited file containing the classes for\n");
    printf("                 the expression data.\n");
    printf(" -c  [<string>]  covariate file. A tab delimited file containing the covariates for\n");
    printf("                 the expression data.\n");
    printf(" -o  [<string>]  output file. If no file is specified it defaults to stdout.\n");
    printf(" -s  [<integer>] number of samples. It should correspond to the number of names\n");
    printf("                 in the first line of the data file.\n");
    printf(" -g  [<integer>] number of variables. It should correspond to the number of lines + 1\n");
    printf("                 in the data file.\n");
    printf(" -n  <string>    variable name. Curve fitting only for the variable name(s)\n");
    printf(" -m  <string>    fit factor name to model (must match one name in the covariate file).\n");
    printf(" -a  <string>    aggregate factor name (must be null or match one name in the factor file).\n");
    printf(" -z  <string>    zero level (must be null or match one level in the aggregate factor).\n");
    printf(" -l  <switch>    log transform the aggregate factor.\n");
    printf(" -t  [string]    type of curve fitting. (default is nlls4p)\n");
    printf(" -e  [<float>]   value that identifies missing data.\n");
    printf(" -p  <switch>    output progress coeeficients.\n\n");

    printf("Examples\n");
    printf("\tfit -d ./t/fit.dat -f ./t/fit.fac -c ./t/fit.cov -m Dose -a Treatment -t nlls4p\n");
    printf("\n");
    exit(0);
  }
   
}

RESULT *
set_result(GROUP *group,
	   int n)
{

  int i, j, idx;
  char str[1000];
  RESULT *res;
  VECTOR *vr;
  char * param[5];

  param[0] = " min";
  param[1] = " max";
  param[2] = " ec50";
  param[3] = " slope";
  param[4] = " pval";
  
  res = TYPE_ALLOC(RESULT);
  res->n = n;
  res->a = group->n * 5;
  res->attributes = ARRAY_ALLOC(res->a, char *);
  res->variables = ARRAY_ALLOC(res->n, char *);
  res->vectors = ARRAY_ALLOC(res->n, VECTOR *);
  for (i = 0; i< n; i++) {
    vr = res->vectors[i] = TYPE_ALLOC(VECTOR);
    vr->data = ARRAY_ALLOC(res->a, double);
    vr->n = res->a;
  }

  idx = 0;
  for (i = 0; i < group->n; i++) {
    for (j = 0; j < 5; j++) {
      strcpy(str, group->names[i]);
      strcat(str, param[j]);
      set_attribute_result(idx,str,res);
      idx++;
    }
  }

  return res;

}

void
set_attribute_result(int i,
		     char *str,
		     RESULT *res)
{

  int len;

  len = strlen(str);
  res->attributes[i] = ARRAY_ALLOC((len+1), char);  
  strncpy(res->attributes[i], str, len);
  res->attributes[i][len] = '\0';

}

RESULT *
fit_data(XDATA *xdata,
	 FACTOR *factor,
	 COVARIATE *covariate,
	 char **variable_names,
	 int variable_number,
	 char *fit_factor,
	 char *agg_factor,
	 char *zero,
	 int logtrans,
	 enum FUNCS func,
	 double *params,
	 double md,
	 int progress)
{

  int i, ii, n, g, gr, len;
  RESULT *res;
  VECTOR *v;
  GROUP *group;

  // Set the number of vectors to stor the data
  if (variable_number > 0) {
    n = variable_number;
  } else {
    n = xdata->n;
  }

  group = set_group(factor, covariate, fit_factor, agg_factor, zero, logtrans);

  res = set_result(group, n);
  
  if (variable_number > 0) {

    ii = 0;
    for (i=0;i<n;i++) {
      g = search_variable(xdata, variable_names[i]);
      if (g < 0) {
	fprintf(stderr, "Unable to find %s in data\n", variable_names[i]);
	exit(1);
      }
      len = strlen(xdata->variables[g]);
      res->variables[ii] = ARRAY_ALLOC((len+1), char);
      strncpy(res->variables[ii], xdata->variables[g], len);
      res->variables[ii][len] = '\0';
      if (verbose || debug) printf("%s\n", xdata->variables[g]);
      v = xdata->vectors[g];
      for (gr = 0; gr < group->n; gr++) {
	if (func == NLLS4P) {
	  fit_nlsr4p(v, group, params, res, ii, gr, progress);
	} else {
	  fit_nlsr4p(v, group, params, res, ii, gr, progress);
	}
      }
      free(v);
      ii++;
    }

  } else {

    for (i=0;i<xdata->n;i++) {
      len = strlen(xdata->variables[i]);
      res->variables[i] = ARRAY_ALLOC((len+1), char);
      strncpy(res->variables[i], xdata->variables[i], len);
      res->variables[i][len] = '\0';
      if (verbose || debug) printf("%s\n", xdata->variables[i]);
      v = xdata->vectors[i];
      for (gr = 0; gr < group->n; gr++) {
	if (func == NLLS4P) {
	  fit_nlsr4p(v, group, params, res, i, gr, progress);
	} else {
	  fit_nlsr4p(v, group, params, res, i, gr, progress);
	}
      }
      free(v);
    }

  }

  free(group);

  return res;

}

void
fit_nlsr4p(VECTOR *v,
	   GROUP *group,
	   double *params,
	   RESULT *res,
	   int idx,
	   int gr,
	   int progress)
{

  int p = 4;

  int i, k, cont, skip, *ia, *sidx;
  double sst, mean, my, dfm, dft, dfe, ssm, msm, mse, f, pval; 
  double *ys, *sig, half, alamda, chisq, ochisq, **covar, **alpha;
  static double *a;
  VECTOR *vr;
  int mrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[],
	      int ma, double **covar, double **alpha, double *chisq,
	      void (*funcs)(double, double [], double *, double [], int), double *alamda);

  a = ARRAY_ALLOC(p, double);
  ia = ARRAY_ALLOC(p, int);

  // Allocate memory
  covar = ARRAY_ALLOC(p, double *);
  alpha = ARRAY_ALLOC(p, double *);
  for (i = 0; i < p; i++) {
    covar[i] = ARRAY_ALLOC(p, double);
    alpha[i] = ARRAY_ALLOC(p, double);    
  }

  // To store results
  vr = res->vectors[idx];

  // Set the fixed parameters
  for (i = 0; i < p; i++) {
    if (params[i] > -999999) {
      a[i] = params[i];
      ia[i] = 0;
    } else {
      a[i] = 0;
      ia[i] = 1;      
    }
  }

  mean = 0.0;
  ys = ARRAY_ALLOC(group->samples[gr], double);
  sig = ARRAY_ALLOC(group->samples[gr], double);
  // Load the x vector by group
  for (i = 0; i < group->samples[gr]; i++) {
    ys[i] = v->data[(group->groups[gr][i])];
    sig[i] = 1.0;
    mean += ys[i];
  }
  mean /= group->samples[gr];
  sst = 0.0;
  for (i = 0; i < group->samples[gr]; i++) {
    my = ys[i] - mean;
    sst += my * my * sig[i];
  }

  // Set the degrees of freedom
  dfm = 4 - 1;
  dft = group->samples[gr] - 1;
  dfe = group->samples[gr] - 4;

  // Set defaults
  skip = 0;
  sidx = ARRAY_ALLOC(group->samples[gr], int);
  sort_by_index(group->samples[gr], ys, sidx);
  if (ia[0] == 1) {
    a[0] = ys[sidx[0]]; 
  }
  if (ia[1] == 1) {
    a[1] = ys[sidx[group->samples[gr]-1]]; 
  }
  if (ia[2] == 1) {
    half = fabs(a[0] + a[1]) / 2;
    for (i = 1; ys[sidx[i]] < half; i++) {
      a[2] = fabs(group->values[gr][sidx[i]] - group->values[gr][sidx[i-1]]) / 2;
    }
    if (group->values[gr][sidx[i]] - group->values[gr][sidx[i-1]] == 0 || a[2] == 0) {
      skip = 1;
    }
  }
  if (ia[3] == 1) {
    a[3] = 1;
  }

  if (verbose || debug) {
    printf("\nLevel %s\n", group->names[gr]);
    printf ("     X => ");
    print_array_d(ys, group->samples[gr]);
    printf ("     Y => ");
    print_array_d(group->values[gr], group->samples[gr]);
    printf ("PARAMS => ");
    print_array_d(a, p);
    printf ("%15s %5s %14s %15s %15s %15s %15s %15s %15s\n",
	    "Variable", "iter", "min", "max", "ec50", "slope", "Chi2", "Lamda", "pval");
  }

  // Check we have more data than parameters

  if (group->levels[gr] > p) {

    // Set the nonlinear fitting

    // Initialization
    alamda = -1;
    if (skip || ! mrqmin(group->values[gr], ys, sig, group->samples[gr],
			 a, ia, p,
			 covar, alpha, &chisq, nlls4p, &alamda)) {
      no_fit(vr, gr, res->variables[idx]);
      free(a);
      free(ia);
      free(covar);
      free(alpha);
      free(ys);
      free(sig);
      free(sidx);
      return;
    }

    // Minimization
    k = 1;
    cont = 1;
    while (k < 50 && cont) {
      ssm = sst - chisq;
      msm = ssm / dfm; 
      mse = chisq / dfe;
      f = fabs(msm/mse);
      pval = betai(0.5*dfe, 0.5*dfm, dfe/(dfe+dfm*f));
      if (verbose || debug) {
	printf("%15s\t%5i\t", res->variables[idx], k);
	printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\n", a[0], a[1], a[2], a[3], chisq, alamda, pval);
      }
      ochisq = chisq;
      mrqmin(group->values[gr], ys, sig, group->samples[gr],
	     a, ia, p,
	     covar, alpha, &chisq, nlls4p, &alamda);
      if ((k > 5) && ((ochisq - chisq) < 0.1)) {
	cont = 0;
      }
      k++;
    }

    // Exit
    alamda = 0.0;
    mrqmin(group->values[gr], ys, sig, group->samples[gr],
	   a, ia, p,
	   covar, alpha, &chisq, nlls4p, &alamda);
    sort_by_index(group->samples[gr], group->values[gr], sidx);
    ssm = sst - chisq;
    msm = ssm / dfm; 
    mse = chisq / dfe;
    f = fabs(msm/mse);
    pval = betai(0.5*dfe, 0.5*dfm, dfe/(dfe+dfm*f)); 
    vr->data[(gr*5)+0] = a[0];
    vr->data[(gr*5)+1] = a[1];
    vr->data[(gr*5)+2] = a[2];
    vr->data[(gr*5)+3] = a[3];
    vr->data[(gr*5)+4] = pval;

  } else {

    no_fit(vr, gr, res->variables[idx]);

  }

  free(a);
  free(ia);
  free(covar);
  free(alpha);
  free(ys);
  free(sig);
  free(sidx);

}

void
no_fit(VECTOR *vr,
       int gr,
       char *variable_name)
{

  vr->data[(gr*5)+0] = MD;
  vr->data[(gr*5)+1] = MD;
  vr->data[(gr*5)+2] = MD;
  vr->data[(gr*5)+3] = MD;
  vr->data[(gr*5)+4] = 1;

  if (verbose || debug) {
    printf("%15s\t%s\n\n", variable_name, "Not applicable");
  }

}

int
search_variable(XDATA *xdata,
		char *variable_name)
{

  int i;
  
  for (i=0;i<xdata->n;i++) {
    if (strcmp(variable_name, xdata->variables[i])==0) {
      return i;
    }
  }
  
  fprintf(stderr, "Unable to find %s in data\n", variable_name);
  exit(1);

}
	
int 
main(int argc,
     char **argv)

{

  char *data_file, *factor_file, *covariate_file, *output_file, **variable_names;
  int  samples, variables, logtrans, progress, lines, variable_number;
  char *fit_factor, *agg_factor, *zero; 
  enum FUNCS func;
  double md, *params;
  XDATA *xdata;
  FACTOR *factor;
  COVARIATE *covariate;
  RESULT *res;
  
  verbose = 0;

  debug = 0;

  parse_arguments(argc, argv,
		  &data_file,
		  &factor_file, &covariate_file,
		  &output_file,
		  &samples, &variables,
		  &variable_names, &variable_number,
		  &fit_factor, &agg_factor, &zero,
		  &logtrans, &func, &params,
		  &md, &progress);
    
  if (samples == 0) {
    samples = count_fields(factor_file);
  }

  lines = count_lines(factor_file);
  lines--;
  
  factor = read_factor(factor_file, lines, samples);

  lines = count_lines(covariate_file);
  lines--;
  
  covariate = read_covariate(covariate_file, lines, samples);

  if (variables == 0) {
    variables = count_lines(data_file);
    variables--;
  }
  
  xdata = read_xdata(data_file, samples, variables);

  res = fit_data(xdata, factor, covariate,
		 variable_names, variable_number,
		 fit_factor, agg_factor, zero, logtrans, func,
		 params, md, progress);

  print_results(output_file, res, 1, 0);
  
  free(data_file);
  free(factor_file);
  free(covariate_file);
  free(output_file);
  free(variable_names);
  free(fit_factor);
  free(agg_factor);
  free(xdata);
  free(factor);
  free(covariate);
  free(res);

  exit(0);

}
