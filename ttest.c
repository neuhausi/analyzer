/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: ttest.c,v 1.9 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#define EXTERN
#include "ttest.h"

void
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **factor_file,
		char **output_file,
		int *samples,
		int *variables,
		enum TT *tt,
		double *md,
		int *permuts)
{

  int c;
  extern char *optarg;
  int errflg = 0;

  *data_file = NULL;

  *factor_file = NULL;

  *output_file = NULL;

  *samples = 0;

  *variables = 0;

  *md = MD;

  *permuts = 10;

  *tt = EQUAL;

  while ((c = getopt(argc, argv, "d:f:o:s:g:t:m:p:vD")) != EOF)
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
    case 't':
      if (strcmp(optarg, "paired") == 0) {
	*tt = PAIRED;
      } else if (strcmp(optarg, "unequal") == 0) {
	*tt = UNEQUAL;
      } else if (strcmp(optarg, "equal") == 0) {
	*tt = EQUAL;
      }
      break;
    case 'm':
      *md = atof(optarg);
      break;
    case 'p':
      *permuts = atoi(optarg);
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
    
  if (*data_file == NULL || *factor_file == NULL ) { 
    errflg++;
  }

  if (errflg) {
    printf("\nUsage: ttest <options> (options in square brackets are optional)\n\n");
    printf(" -d  <string>    data file. A tab delimeted file with the following format.\n");
    printf("                 A header line with sample names. Lines with expression values\n");
    printf("                 preceded by the variable name at the begining of the line.\n");
    printf(" -f  <string>    factor file. A tab delimet file containing the classes for\n");
    printf("                 the expression data.\n");
    printf(" -o  [<string>]  output file. If no file is specified it defaults to stdout.\n");
    printf(" -s  [<integer>] number of samples. It should correspond to the number of names\n");
    printf("                 in the first line of the data file.\n");
    printf(" -g  [<integer>] number of variables. It should correspond to the number of lines + 1\n");
    printf("                 in the data file.\n");
    printf(" -t  [<string>]  [equal|unequal|paired]. Make reference to the variance in the\n");
    printf("                 groups or whether to use a paired ttest.\n");
    printf("                 (default is 'equal' variance).\n");
    printf(" -m  [<float>]   value that identifies missing data.\n");
    printf(" -p  [<integer>] number of random permutations (not implemented yet).\n\n");
    printf("Examples\n");
    printf("\n\t# t-test (equal variance)\n");
    printf("\tttest -d ./t/ttest.dat -f ./t/ttest.fac -t equal\n");
    printf("\n\t# t-test (unequal variance)\n");
    printf("\tttest -d ./t/ttest.dat -f ./t/ttest.fac -t unequal\n");
    printf("\n\t# t-test (paired)\n");
    printf("\tttest -d ./t/ttest.dat -f ./t/ttest.fac -t paired\n");
    printf("\n");
    exit(0);
  }
   
}

RESULT *
run_ttest(XDATA *xdata,
	  FACTOR *factor,
	  enum TT tt, 
	  double md,
	  int permuts)
{

  int i, j, len, *m_cnt;
  double p, t, df, m1, m2, v1, v2, *data, *data1, *data2;
  char str[10000];
  VECTOR *v, *vr;
  RESULT *res;

  m_cnt = ARRAY_ALLOC(factor->classes[0], int);
  data1 = ARRAY_ALLOC(factor->members[0][0], double);
  data2 = ARRAY_ALLOC(factor->members[0][1], double);

  res = TYPE_ALLOC(RESULT);
  res->n = xdata->n;
  res->a = 5;
  res->attributes = ARRAY_ALLOC(res->a, char *);
  res->variables = ARRAY_ALLOC(res->n, char *);
  res->vectors = ARRAY_ALLOC(res->n, VECTOR *);
  strcpy(str, factor->fnames[0]);
  len = strlen(str);
  res->attributes[0] = ARRAY_ALLOC((len+1), char);
  strncpy(res->attributes[0], str, len);
  res->attributes[0][len] = '\0';
  strcpy(str, "Mean ");
  strcat(str, factor->cnames[0][0]);
  len = strlen(str);
  res->attributes[1] = ARRAY_ALLOC((len+1), char);
  strncpy(res->attributes[1], str, len);
  res->attributes[1][len] = '\0';
  strcpy(str, "Var ");
  strcat(str, factor->cnames[0][0]);
  len = strlen(str);
  res->attributes[2] = ARRAY_ALLOC((len+1), char);
  strncpy(res->attributes[2], str, len);
  res->attributes[2][len] = '\0';
  strcpy(str, "Mean ");
  strcat(str, factor->cnames[0][1]);
  len = strlen(str);
  res->attributes[3] = ARRAY_ALLOC((len+1), char);
  strncpy(res->attributes[3], str, len);
  res->attributes[3][len] = '\0';
  strcpy(str, "Var ");
  strcat(str, factor->cnames[0][1]);
  len = strlen(str);
  res->attributes[4] = ARRAY_ALLOC((len+1), char);
  strncpy(res->attributes[4], str, len);
  res->attributes[4][len] = '\0';

  m_cnt[0] = 0;
  m_cnt[1] = 0;


  if (md > MD) {

    for (i = 0; i < xdata->n; i++) {

      len = strlen(xdata->variables[i]);
      res->variables[i] = ARRAY_ALLOC((len+1), char);
      strncpy(res->variables[i], xdata->variables[i], len);
      res->variables[i][len] = '\0';

      vr = res->vectors[i] = TYPE_ALLOC(VECTOR);
      vr->data = ARRAY_ALLOC(res->a, double);
      vr->n = res->a;

      v = xdata->vectors[i];
      data = v->data;
      m_cnt[0] = 0;
      m_cnt[1] = 0;
      for (j = 0; j < xdata->s; j++) {
	if (factor->factors[0][j] == 0 && data[j] != md) {
	  data1[m_cnt[0]] = data[j];
	  m_cnt[0]++;
	} else if (factor->factors[0][j] == 1 && data[j] != md) {
	  data2[m_cnt[1]] = data[j];
	  m_cnt[1]++;
	}
      }
      if (tt == EQUAL) {
	ttest(data1, m_cnt[0], data2, m_cnt[1], &t, &df, &p,
	      &m1, &v1, &m2 , &v2);
      } else if (tt == UNEQUAL) {
	tutest(data1, m_cnt[0], data2, m_cnt[1], &t, &df, &p,
	       &m1, &v1, &m2 , &v2);
      } else if (tt == PAIRED) {
	if (m_cnt[0] != m_cnt[1]) {
	  exit_on_error("Not the same nuumber of values");
	}
	tptest(data1, data2, m_cnt[0], &t, &df, &p,
	       &m1, &v1, &m2 , &v2);
      }
      vr->data[0] = p;
      vr->data[1] = m1;
      vr->data[2] = v1;
      vr->data[3] = m2;
      vr->data[4] = v2;
    }

  } else {

    for (i = 0; i < xdata->n; i++) {

      len = strlen(xdata->variables[i]);
      res->variables[i] = ARRAY_ALLOC((len+1), char);
      strncpy(res->variables[i], xdata->variables[i], len);
      res->variables[i][len] = '\0';

      vr = res->vectors[i] = TYPE_ALLOC(VECTOR);
      vr->data = ARRAY_ALLOC(res->a, double);
      vr->n = res->a;

      v = xdata->vectors[i];
      data = v->data;
      m_cnt[0] = 0;
      m_cnt[1] = 0;
      for (j = 0; j < xdata->s; j++) {
	if (factor->factors[0][j] == 0) {
	  data1[m_cnt[0]] = data[j];
	  m_cnt[0]++;
	} else if (factor->factors[0][j] == 1) {
	  data2[m_cnt[1]] = data[j];
	  m_cnt[1]++;
	}
      }
      if (tt == EQUAL) {
	ttest(data1, factor->members[0][0], data2, factor->members[0][1],
	      &t, &df, &p, &m1, &v1, &m2 , &v2);
      } else if (tt == UNEQUAL) {
	tutest(data1, factor->members[0][0], data2, factor->members[0][1],
	       &t, &df, &p, &m1, &v1, &m2 , &v2);
      } else if (tt == PAIRED) {
	if (factor->members[0][0] != factor->members[0][1]) {
	  exit_on_error("Not the same nuumber of values");
	}
	tptest(data1, data2, factor->members[0][0],
	       &t, &df, &p, &m1, &v1, &m2 , &v2);
      }
      vr->data[0] = p;
      vr->data[1] = m1;
      vr->data[2] = v1;
      vr->data[3] = m2;
      vr->data[4] = v2;
    }

  }

  return res;

}

void
ttest(double data1[], int n1, double data2[], int n2, double *t, double *df,
      double *prob, double *m1, double *v1, double *m2, double *v2)
{

  double var1,var2,svar,ave1,ave2;
  
  avevar(data1,n1,&ave1,&var1);
  avevar(data2,n2,&ave2,&var2);
  if (debug) {
      printf ("Data1 (%i) => ", n1);
      print_array_d(data1, n1);
      printf ("\n");
      printf ("Data2 (%i) => ", n2);
      print_array_d(data2, n2);
      printf ("\n");
  }
  *m1 = ave1;
  *v1 = var1;
  *m2 = ave2;
  *v2 = var2;
  *df=n1+n2-2;
  svar=((n1-1)*var1+(n2-1)*var2)/(*df);
  *t=(ave1-ave2)/sqrt(svar*(1.0/n1+1.0/n2));
  *prob=betai(0.5*(*df),0.5,(*df)/((*df)+(*t)*(*t)));

}

void
tutest(double data1[], int n1, double data2[], int n2, double *t, double *df,
       double *prob, double *m1, double *v1, double *m2, double *v2)
{

  double var1,var2,ave1,ave2;
  
  avevar(data1,n1,&ave1,&var1);
  avevar(data2,n2,&ave2,&var2);
  if (debug) {
      printf ("Data1 (%i) => ", n1);
      print_array_d(data1, n1);
      printf ("\n");
      printf ("Data2 (%i) => ", n2);
      print_array_d(data2, n2);
      printf ("\n");
  }
  *m1 = ave1;
  *v1 = var1;
  *m2 = ave2;
  *v2 = var2;
  *t=(ave1-ave2)/sqrt(var1/n1+var2/n2);
  *df=SQR(var1/n1+var2/n2)/(SQR(var1/n1)/(n1-1)+SQR(var2/n2)/(n2-1));
  *prob=betai(0.5*(*df),0.5,(*df)/((*df)+SQR(*t)));

}

void
tptest(double data1[], double data2[], int n, double *t, double *df,
       double *prob, double *m1, double *v1, double *m2, double *v2)
{

  int j;
  double var1,var2,ave1,ave2,sd,cov=0.0;
  
  avevar(data1,n,&ave1,&var1);
  avevar(data2,n,&ave2,&var2);
  if (debug) {
      printf ("Data1 (%i) => ", n);
      print_array_d(data1, n);
      printf ("\n");
      printf ("Data2 (%i) => ", n);
      print_array_d(data2, n);
      printf ("\n");
  }
  *m1 = ave1;
  *v1 = var1;
  *m2 = ave2;
  *v2 = var2;
  for (j=0;j<n;j++)
    cov += (data1[j]-ave1)*(data2[j]-ave2);
  cov /= *df=n-1;
  sd=sqrt((var1+var2-2.0*cov)/n);
  *t=(ave1-ave2)/sd;
  *prob=betai(0.5*(*df),0.5,(*df)/((*df)+(*t)*(*t)));

}

int
main(int argc,
     char **argv)

{

  char *data_file, *factor_file, *output_file;
  int  samples, variables, permuts, lines;
  double md;
  enum TT tt;
  XDATA *xdata;
  FACTOR *factor;
  RESULT *res;

  verbose = 0;

  debug = 0;

  parse_arguments(argc, argv, &data_file, &factor_file, &output_file,
		  &samples, &variables, &tt, &md, &permuts);

  if (samples == 0) {
    samples = count_fields(factor_file);
  }

  lines = 1;

  factor = read_factor(factor_file, lines, samples);

  if (variables == 0) {
    variables = count_lines(data_file);
    variables--;
  }

  xdata = read_xdata(data_file, samples, variables);

  res = run_ttest(xdata, factor, tt, md, permuts);

  print_results(output_file, res, 1, 0);

  exit(0);

}
