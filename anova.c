/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: anova.c,v 1.9 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#define EXTERN
#include "anova.h"

void
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **factor_file,
		char **output_file,
		int *samples,
		int *variables,
		int *inter,
		enum TT *tt,
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
  
  *inter = 1;

  *tt = BETWEEN;
  
  *permuts = 10;
  
  while ((c = getopt(argc, argv, "d:f:o:s:g:nt:p:vD")) != EOF)
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
    case 'n':
      *inter = 0;
      break;
    case 't':
      if (strcmp(optarg, "within") == 0) {
	*tt = WITHIN;
      } else if (strcmp(optarg, "mixed") == 0) {
	*tt = MIXED;
      } else {
	*tt = BETWEEN;
      }
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
    printf("\nUsage: anova <options>  (options in square brackets are optional)\n\n");
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
    printf(" -t  [<string>]  [between|within|mixed]. Type of ANOVA. For within and mixed\n");
    printf("                 models the subject factor must be the first factor in the file.\n");
    printf("                 (default is between).\n");
    printf(" -n  [switch]    No interaction. Defaults to include interaction between factors.\n");
    printf(" -p  [<integer>] number of random permutations (not implemented yet).\n\n");
    printf("Examples\n");
    printf("\n\t# One-way ANOVA (balanced)\n");
    printf("\tanova -d ./t/1f.3l.b.dat -f ./t/1f.3l.b.fac\n");
    printf("\n\t# Three-way ANOVA (missing data no interaction)\n");
    printf("\tanova -d ./t/3f.2l.2l.2l.md.dat -f ./t/3f.2l.2l.2l.b.fac -n\n");
    printf("\n\t# One factor within subject ANOVA\n");
    printf("\tanova -d ./t/2f.6l.2l.w.dat -f ./t/2f.6l.2l.w.fac -t within\n");
    printf("\n\t# Two factor mixed within subject ANOVA\n");
    printf("\tanova -d ./t/3f.4l.3l.4l.m.dat -f ./t/3f.4l.3l.4l.m.fac -t mixed\n");
    printf("\n");
    exit(0);
  }
   
}

RESULT *
run_one_way_anova(XDATA *xdata,
		  FACTOR *factor, 
		  int permuts)
{

  int i, j, c, cnt, tot, len, *cfactor;
  double p, *data, *cdata, *cls_avg;
  char str[10000];
  VECTOR *v, *vr;
  RESULT *res;

  c = factor->classes[0];
  res = TYPE_ALLOC(RESULT);
  res->n = xdata->n;
  res->a = 1 + c;
  res->attributes = ARRAY_ALLOC(res->a, char *);
  res->variables = ARRAY_ALLOC(res->n, char *);
  res->vectors = ARRAY_ALLOC(res->n, VECTOR *);
  
  for (i = 0; i < c; i++) {
    strcpy(str, "Mean ");
    strcat(str, factor->cnames[0][i]);
    len = strlen(str);
    res->attributes[i] = ARRAY_ALLOC((len+1), char);
    strncpy(res->attributes[i], str, len);
    res->attributes[i][len] = '\0';
  }
  strcpy(str, factor->fnames[0]);
  len = strlen(str);
  res->attributes[c] = ARRAY_ALLOC((len+1), char);
  strncpy(res->attributes[c], str, len);
  res->attributes[c][len] = '\0';

  tot = 0;
  cnt = 0;
  for (i = 0; i < c; i++) {
    tot += factor->members[0][i];
  }
  cdata   = ARRAY_ALLOC(tot, double);
  cfactor = ARRAY_ALLOC(tot, int);
  cls_avg = ARRAY_ALLOC(c, double);

  for (i = 0; i < xdata->n; i++) {

    len = strlen(xdata->variables[i]);
    res->variables[i] = ARRAY_ALLOC((len+1), char);
    strncpy(res->variables[i], xdata->variables[i], len);
    res->variables[i][len] = '\0';

    vr = res->vectors[i] = TYPE_ALLOC(VECTOR);
    vr->data = ARRAY_ALLOC(res->a, double);
    vr->n = res->a;

    cnt = 0;
    v = xdata->vectors[i];
    data = v->data;
    for (j = 0; j < xdata->s; j++) {
      if (factor->factors[0][j] >= 0 && ! isnan(data[j])) {
	cdata[cnt] = data[j];
	cfactor[cnt] = factor->factors[0][j];
	cnt++;
      }
    }
    if (verbose || debug) printf("%s\n", xdata->variables[i]);
    anova_one_way(cdata, cfactor, cnt, c, &p, &factor->fnames[0], &(*cls_avg));
    for (j = 0; j < c; j++) {
      vr->data[j] = cls_avg[j];
    }
    vr->data[c] = p;
  }

  return res;

}

RESULT *
run_multi_factorial_anova(XDATA *xdata,
			  FACTOR *factor, 
			  enum TT tt, 
			  int inter,
			  int permuts)
{

  int i, j, nf, cnt, tot, len, *cls, **cfactor;
  double p1, p2, p3, p12, p13, p23, p123, *data, *cdata;
  char str[1000];
  VECTOR *v, *vr;
  RESULT *res;

  nf = factor->n;
  cls = ARRAY_ALLOC(nf, int);
  cfactor = ARRAY_ALLOC(nf, int *);
  res = TYPE_ALLOC(RESULT);
  res->n = xdata->n;
  res->a = nf == 2 ? 3 : 7;
  if (tt == BETWEEN) {
    if (! inter) {
      res->a = nf == 2 ? 2 : 3;
    }
    res->attributes = ARRAY_ALLOC(res->a, char *);
    res->variables = ARRAY_ALLOC(res->n, char *);
    res->vectors = ARRAY_ALLOC(res->n, VECTOR *);
    strcpy(str, factor->fnames[0]);
    len = strlen(str);
    res->attributes[0] = ARRAY_ALLOC((len+1), char);
    strncpy(res->attributes[0], str, len);
    res->attributes[0][len] = '\0';
    strcpy(str, factor->fnames[1]);
    len = strlen(str);
    res->attributes[1] = ARRAY_ALLOC((len+1), char);
    strncpy(res->attributes[1], str, len);
    res->attributes[1][len] = '\0';
    if (nf==2 && inter) {
      strcpy(str, factor->fnames[0]);
      strcat(str, " * ");
      strcat(str, factor->fnames[1]);
      len = strlen(str);
      res->attributes[2] = ARRAY_ALLOC((len+1), char);
      strncpy(res->attributes[2], str, len);
      res->attributes[2][len] = '\0';
    } else if (nf==3) {
      strcpy(str, factor->fnames[2]);
      len = strlen(str);
      res->attributes[2] = ARRAY_ALLOC((len+1), char);
      strncpy(res->attributes[2], str, len);
      res->attributes[2][len] = '\0';
      if (inter) {
	strcpy(str, factor->fnames[0]);
	strcat(str, " * ");
	strcat(str, factor->fnames[1]);
	len = strlen(str);
	res->attributes[3] = ARRAY_ALLOC((len+1), char);
	strncpy(res->attributes[3], str, len);
	res->attributes[3][len] = '\0';
	strcpy(str, factor->fnames[0]);
	strcat(str, " * ");
	strcat(str, factor->fnames[2]);
	len = strlen(str);
	res->attributes[4] = ARRAY_ALLOC((len+1), char);
	strncpy(res->attributes[4], str, len);
	res->attributes[4][len] = '\0';
	strcpy(str, factor->fnames[1]);
	strcat(str, " * ");
	strcat(str, factor->fnames[2]);
	len = strlen(str);
	res->attributes[5] = ARRAY_ALLOC((len+1), char);
	strncpy(res->attributes[5], str, len);
	res->attributes[5][len] = '\0';
	strcpy(str, factor->fnames[0]);
	strcat(str, " * ");
	strcat(str, factor->fnames[1]);
	strcat(str, " * ");
	strcat(str, factor->fnames[2]);
	len = strlen(str);
	res->attributes[6] = ARRAY_ALLOC((len+1), char);
	strncpy(res->attributes[6], str, len);
	res->attributes[6][len] = '\0';
      }
    }      
  } else if (tt == WITHIN || tt == MIXED) {
    res->a = nf == 2 ? 1 : 3;
    res->attributes = ARRAY_ALLOC(res->a, char *);
    res->variables = ARRAY_ALLOC(res->n, char *);
    res->vectors = ARRAY_ALLOC(res->n, VECTOR *);
    strcpy(str, factor->fnames[1]);
    len = strlen(str);
    res->attributes[0] = ARRAY_ALLOC((len+1), char);
    strncpy(res->attributes[0], str, len);
    res->attributes[0][len] = '\0';
    if (nf==3) {
      strcpy(str, factor->fnames[2]);
      len = strlen(str);
      res->attributes[1] = ARRAY_ALLOC((len+1), char);
      strncpy(res->attributes[1], str, len);
      res->attributes[1][len] = '\0';
      strcpy(str, factor->fnames[1]);
      strcat(str, " * ");
      strcat(str, factor->fnames[2]);
      len = strlen(str);
      res->attributes[2] = ARRAY_ALLOC((len+1), char);
      strncpy(res->attributes[2], str, len);
      res->attributes[2][len] = '\0';
    }
  }

  for (i = 0; i < factor->n; i++) {
    tot = 0;
    cls[i] = factor->classes[i];
    for (j = 0; j < factor->classes[i]; j++) {
      tot += factor->members[i][j];
    }
    cfactor[i] = ARRAY_ALLOC(tot, int);
  }
  cdata   = ARRAY_ALLOC(tot, double);
  
  for (i = 0; i < xdata->n; i++) {

    len = strlen(xdata->variables[i]);
    res->variables[i] = ARRAY_ALLOC((len+1), char);
    strncpy(res->variables[i], xdata->variables[i], len);
    res->variables[i][len] = '\0';

    vr = res->vectors[i] = TYPE_ALLOC(VECTOR);
    vr->data = ARRAY_ALLOC(res->a, double);
    vr->n = res->a;

    cnt = 0;
    v = xdata->vectors[i];
    data = v->data;

    for (j = 0; j < xdata->s; j++) {
      if (nf==2 &&
	factor->factors[0][j] >= 0 &&
	factor->factors[1][j] >= 0 &&
	! isnan(data[j])) {
	if (factor->factors[0][j] >= 0 && factor->factors[1][j] >= 0 ) {
	  cdata[cnt] = data[j];
	  cfactor[0][cnt] = factor->factors[0][j];
	  cfactor[1][cnt] = factor->factors[1][j];
	  cnt++;
	}
      } else if (nf==3 &&
		 factor->factors[0][j] >= 0 &&
		 factor->factors[1][j] >= 0 &&
		 factor->factors[2][j] >= 0 &&
		 ! isnan(data[j])) {
	cdata[cnt] = data[j];
	cfactor[0][cnt] = factor->factors[0][j];
	cfactor[1][cnt] = factor->factors[1][j];
	cfactor[2][cnt] = factor->factors[2][j];
	cnt++;
      }
    }

    if (tt == BETWEEN) {
      if (verbose || debug) printf("%s\n", xdata->variables[i]);
      if (nf==2) {
	anova_two_way(cdata, cfactor, cnt, cls, inter,
		      &p1, &p2, &p12,
		      &factor->fnames[0],  &factor->fnames[1]);
	vr->data[0] = p1;
	vr->data[1] = p2;
	vr->data[2] = p12;
      } else {
	anova_three_way(cdata, cfactor, cnt, cls, inter,
			&p1, &p2, &p3, &p12, &p13, &p23, &p123,
			&factor->fnames[0],  &factor->fnames[1],
			&factor->fnames[2]);
	vr->data[0] = p1;
	vr->data[1] = p2;
	vr->data[2] = p3;
	vr->data[3] = p12;
	vr->data[4] = p13;
	vr->data[5] = p23;
	vr->data[6] = p123;
      }
    } else if (tt == WITHIN || tt == MIXED) {
      if (verbose || debug) printf("%s\n", xdata->variables[i]);
      if (nf==2) {
	anova_one_way_within(cdata, cfactor, cnt, cls,
			     &p1, &factor->fnames[0],
			     &factor->fnames[1]);
	vr->data[0] = p1;
      } else {
	if (tt == MIXED) {
	  anova_two_way_mixed(cdata, cfactor, cnt, cls,
			      &p2, &p3, &p23,
			      &factor->fnames[0],
			      &factor->fnames[1],
			      &factor->fnames[2]);
	  vr->data[0] = p2;
	  vr->data[1] = p3;
	  vr->data[2] = p23;
	} else {
	  anova_two_way_within(cdata, cfactor, cnt, cls,
			       &p2, &p3, &p23,
			       &factor->fnames[0],
			       &factor->fnames[1],
			       &factor->fnames[2]);
	  vr->data[0] = p2;
	  vr->data[1] = p3;
	  vr->data[2] = p23;
	}
      }	  
    }

  }

  return res;

}

void
anova_one_way(double data[], int factor[], int n_tot, int cls,
	      double *prob, char **name, double *cls_avg)
{
    
  // This code utilizes the traditional unweighted means
  
  int i;
  double T, A, Y, n;
  double SST, SSA, SSY;
  double MSA, MSY;
  int dfT, dfA, dfY;
  double f;
  double cell_sum[cls];
  double cell_ssq[cls];
  double cell_avg[cls];
  double cell_uav[cls];
  double cell_mem[cls];
  char	str[25];
  
  // Initialize the arrays otherwise weird things may happen
  
  
  for (i=0;i<cls;i++) {
    cell_sum[i] = 0;
    cell_ssq[i] = 0;
    cell_avg[i] = 0;
    cell_uav[i] = 0;
    cell_mem[i] = 0;
    cls_avg[i]  = 0;
  }
  
  // Calculate cell sums and cell sum of squares
  
  for (i=0;i<n_tot;i++) {
    cell_ssq[factor[i]] += data[i]*data[i]; 
    cell_sum[factor[i]] += data[i]; 
    cell_mem[factor[i]]++;
  }
  
  // Calculate cell averages, harmonic mean, within groups degrees of freedom
  // and sum of squares for unbalanced designs
  
  n=0;
  Y=0;
  dfY=0;
  for (i=0;i<cls;i++) {
    cell_avg[i]=cell_sum[i]/cell_mem[i];
    cls_avg[i]=cell_avg[i];
    n += (1/cell_mem[i]);
    Y += cell_ssq[i]-(cell_sum[i]*cell_sum[i]/cell_mem[i]);
    dfY += (cell_mem[i]-1);
  }
  n=cls/n;
  
  // Set the unweigthed cell averages and calculate the grand sum
  // and factor sum of squares
  
  T=0;
  A=0;
  for (i=0;i<cls;i++) {
    cell_uav[i]=cell_avg[i]*n;
    T += cell_uav[i];
    A += cell_uav[i]*cell_uav[i];
  }
  T=(T*T)/(cls*n);
  A /= n;
  
  // Sum of Squares
  SSA = A-T+TINY;
  SSY = Y+TINY;
  SST = SSA+SSY+TINY;
  
  // Degrees of freedom
  
  dfA = cls-1+TINY;
  dfT = n_tot-1;
  
  // Mean sum of Squares
  
  MSA=SSA/dfA;
  MSY=SSY/dfY;
  
  // F-ratios
  
  f=MSA/MSY;
  
  // Probability
  
  *prob=betai(0.5*dfY,0.5*dfA,dfY/(dfY+dfA*(f)));
  
  if (verbose || debug) {
    printf ("SOURCE                               SS     df           MS         F    Pr(>F)\n");
    strcpy(str, *name);
    strcat(str, " (A)");
    printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSA, dfA, MSA, f, *prob);
    strcpy(str, "Within (error)");
    printf ("%-24s %14.2f %6i %12.2f\n", str, SSY, dfY, MSY);
    strcpy(str, "Total");
    printf ("%-24s %14.2f %6i\n\n", str, SST, dfT);
  }
  
}

void
anova_one_way_within(double data[], int **factor, int n_tot, int cls[], double *prob,
		     char **nameS, char **name)
{

  // This code can only support balanced within-subjects designs
  
  int i, j;
  double T, A, S, Y;
  double SST, SSA, SSS, SSY;
  double MSA, MSS, MSY;
  int dfT, dfA, dfS, dfY;
  double fA;
  double col_sum[cls[0]];
  double row_sum[cls[1]];
  double cell_sum[cls[0]][cls[1]];
  double cell_ssq[cls[0]][cls[1]];
  char	str[25];
  
  // Initialize the arrays otherwise weird things may happen
  
  for (i=0;i<cls[0];i++) {
    col_sum[i] = 0;
  }
  for (i=0;i<cls[1];i++) {
    row_sum[i] = 0;    
  }
  for (i=0;i<cls[0];i++) {
    for (j=0;j<cls[1];j++) {
      cell_sum[i][j] = 0;
      cell_ssq[i][j] = 0;
    }
  }
  
  // Calculate cell sums and cell sum of squares
  
  for (i=0;i<n_tot;i++) {
    cell_ssq[factor[0][i]][factor[1][i]] += data[i]*data[i]; 
    cell_sum[factor[0][i]][factor[1][i]] += data[i]; 
  }
  
  // Calculate the grand sum, rows and columns averages
  
  T=0;
  Y=0;
  for (i=0;i<cls[0];i++) {
    for (j=0;j<cls[1];j++) {
      T += cell_sum[i][j];
      Y += cell_ssq[i][j];
      col_sum[i] += cell_sum[i][j];
      row_sum[j] += cell_sum[i][j];
    }
  }
  T=(T*T)/(cls[0]*cls[1]);
  
  // Calculate the factors sum of squares
  
  S=0;
  for (i=0;i<cls[0];i++) {
    S += col_sum[i]*col_sum[i];
  }
  S /= cls[1];
  A=0;
  for (i=0;i<cls[1];i++) {
    A += row_sum[i]*row_sum[i];
  }    
  A /= cls[0];
  
  // Sum of Squares
  
  SSA = A-T+TINY;
  SSS = S-T+TINY;
  SSY = Y-A-S+T+TINY;
  SST = SSA+SSS+SSY+TINY;
  
  // Degrees of freedom
  
  dfA = cls[1]-1+TINY;
  dfS = cls[0]-1+TINY;
  dfY = dfA*dfS;
  dfT = n_tot-1;
  
  // Mean sum of Squares
  
  MSA=SSA/dfA;
  MSS=SSS/dfS;
  MSY=SSY/dfY;
  
  // F-ratios
  
  fA=MSA/MSY;
  
  // Probabilities
  
  *prob=betai(0.5*dfY,0.5*dfA,dfY/(dfY+dfA*(fA)));

  if (verbose || debug) {
    printf ("DATA\n");
    print_anova_data(data, n_tot, factor, 2);
    printf ("\nSOURCE                               SS     df           MS         F    Pr(>F)\n");
    strcpy(str, *name);
    strcat(str, " (A)");
    printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSA, dfA, MSA, fA, *prob);
    strcpy(str, *nameS);
    strcat(str, " (S)");
    printf ("%-24s %14.2f %6i %12.2f\n", str, SSS, dfS, MSS);
    strcpy(str, "A*S");
    printf ("%-24s %14.2f %6i %12.2f\n", str, SSY, dfY, MSY);
    strcpy(str, "Total");
    printf ("%-24s %14.2f %6i\n", str, SST, dfT);
  }
  
}

void
anova_two_way(double data[], int **factor, int n_tot, int cls[],
	      int inter, double *probA, double *probB, double *probAB,
	      char **nameA, char **nameB)
{

  // This code utilizes the traditional unweighted means
  
  int i, j;
  double T, A, B, AB, Y, n;
  double SST, SSA, SSB, SSAB, SSY;
  double MSA, MSB, MSAB, MSY;
  int dfT, dfA, dfB, dfAB, dfY;
  double fA, fB, fAB;
  double col_sum[cls[0]];
  double row_sum[cls[1]];
  double cell_sum[cls[0]][cls[1]];
  double cell_ssq[cls[0]][cls[1]];
  double cell_avg[cls[0]][cls[1]];
  double cell_uav[cls[0]][cls[1]];
  double cell_mem[cls[0]][cls[1]];
  char	str[25];
  
  // Initialize the arrays otherwise weird things may happen
  
  for (i=0;i<cls[0];i++) {
    col_sum[i] = 0;
  }
  for (i=0;i<cls[1];i++) {
    row_sum[i] = 0;    
  }
  for (i=0;i<cls[0];i++) {
    for (j=0;j<cls[1];j++) {
      cell_sum[i][j] = 0;
      cell_ssq[i][j] = 0;
      cell_avg[i][j] = 0;
      cell_uav[i][j] = 0;
      cell_mem[i][j] = 0;
    }
  }
  
  // Calculate cell sums and cell sum of squares
  
  for (i=0;i<n_tot;i++) {
    cell_ssq[factor[0][i]][factor[1][i]] += data[i]*data[i]; 
    cell_sum[factor[0][i]][factor[1][i]] += data[i]; 
    cell_mem[factor[0][i]][factor[1][i]]++;
  }
  
  // Calculate cell averages, harmonic mean, within groups degrees of freedom
  // and sum of squares for unbalanced designs
  
  n=0;
  Y=0;
  dfY=0;
  for (i=0;i<cls[0];i++) {
    for (j=0;j<cls[1];j++) {
      cell_avg[i][j]=cell_sum[i][j]/cell_mem[i][j];
      n += (1/cell_mem[i][j]);
      Y += cell_ssq[i][j]-(cell_sum[i][j]*cell_sum[i][j]/cell_mem[i][j]);
      dfY += (cell_mem[i][j]-1);
    }
  }
  n=cls[0]*cls[1]/n;
  
  // Set the unweigthed cell averages and calculate the grand sum,
  // rows and columns averages, and factor interaction sum of squares
  
  T=0;
  AB=0;
  for (i=0;i<cls[0];i++) {
    for (j=0;j<cls[1];j++) {
      cell_uav[i][j]=cell_avg[i][j]*n;
      T += cell_uav[i][j];
      AB += cell_uav[i][j]*cell_uav[i][j];
      col_sum[i] += cell_uav[i][j];
      row_sum[j] += cell_uav[i][j];
    }
  }
  T=(T*T)/(cls[0]*cls[1]*n);
  AB /= n;
  
  // Calculate the factors sum of squares
  
  A=0;
  for (i=0;i<cls[0];i++) {
    A += col_sum[i]*col_sum[i];
  }
  B=0;
  for (i=0;i<cls[1];i++) {
    B += row_sum[i]*row_sum[i];
  }    
  A /= cls[1]*n;
  B /= cls[0]*n;
  
  // Sum of Squares
  
  SSA = A-T+TINY;
  SSB = B-T+TINY;
  SSAB = AB-A-B+T+TINY;
  if (! inter) {
    SSY = Y+SSAB+TINY;
    SST = SSA+SSB+SSY+TINY;
  } else {
    SSY = Y+TINY;
    SST = SSA+SSB+SSAB+SSY+TINY;
  }
  
  // Degrees of freedom
  
  dfA = cls[0]-1+TINY;
  dfB = cls[1]-1+TINY;
  dfAB = dfA*dfB;
  dfT = n_tot-1;
  if (! inter) dfY=(n_tot-dfA-dfB-1)+TINY;
  
  // Mean sum of Squares
  
  MSA=SSA/dfA;
  MSB=SSB/dfB;
  MSAB=SSAB/dfAB;
  MSY=SSY/dfY;
  
  // F-ratios
  
  fA=MSA/MSY;
  fB=MSB/MSY;
  fAB=MSAB/MSY;
  
  // Probabilities
  
  *probA=betai(0.5*dfY,0.5*dfA,dfY/(dfY+dfA*(fA)));
  *probB=betai(0.5*dfY,0.5*dfB,dfY/(dfY+dfB*(fB)));
  *probAB=betai(0.5*dfY,0.5*dfAB,dfY/(dfY+dfAB*(fAB)));
  
  if (verbose || debug) {
    printf ("DATA\n");
    print_anova_data(data, n_tot, factor, 2);
    printf ("\nSOURCE                               SS     df           MS         F    Pr(>F)\n");
    strcpy(str, *nameA);
    strcat(str, " (A)");
    printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSA, dfA, MSA, fA, *probA);
    strcpy(str, *nameB);
    strcat(str, " (B)");
    printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSB, dfB, MSB, fB, *probB);
    if (inter) {
      strcpy(str, "A*B");
      printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSAB, dfAB, MSAB, fAB, *probAB);
    } else {
      *probAB=1.0;
    }
    strcpy(str, "Within (error)");
    printf ("%-24s %14.2f %6i %12.2f\n", str, SSY, dfY, MSY);
    strcpy(str, "Total");
    printf ("%-24s %14.2f %6i\n", str, SST, dfT);
  }
  
}

void
anova_two_way_mixed(double data[], int **factor, int n, int cls[],
		    double *probA, double *probB, double *probAB,
		    char **nameS, char **nameA, char **nameB)
{
  
  // This code can only support balanced mixed model designs

  int i, j, k, idx;
  double T, A, B, AB, AS, Y;
  double SST, SSA, SSB, SSAB, SSAS, SSY;
  double MSA, MSB, MSAB, MSAS, MSY;
  int dfT, dfA, dfB, dfAB, dfAS, dfY;
  double fA, fB, fAB;
  double col_sumA[cls[1]];
  double col_sumB[cls[2]];
  double col_row_sumAS[cls[0]*cls[1]];
  double col_row_sumAB[cls[1]*cls[2]];
  double cell_sum[cls[0]][cls[1]][cls[2]];
  double cell_ssq[cls[0]][cls[1]][cls[2]];
  char	str[25];
  
  // Initialize the arrays otherwise weird things may happen
  
  for (i=0;i<cls[1];i++) {
    col_sumA[i] = 0;    
  }
  for (i=0;i<cls[2];i++) {
    col_sumB[i] = 0;    
  }
  for (i=0;i<(cls[0]*cls[1]);i++) {
    col_row_sumAS[i] = 0;
  }
  for (i=0;i<(cls[1]*cls[2]);i++) {
    col_row_sumAB[i] = 0;
  }
  for (i=0;i<cls[0];i++) {
    for (j=0;j<cls[1];j++) {
      for (k=0;k<cls[2];k++) {
	cell_sum[i][j][k] = 0;
	cell_ssq[i][j][k] = 0;
      }
    }
  }
  
  // Calculate cell sums and cell sum of squares
  
  for (i=0;i<n;i++) {
    cell_ssq[factor[0][i]][factor[1][i]][factor[2][i]] += data[i]*data[i]; 
    cell_sum[factor[0][i]][factor[1][i]][factor[2][i]] += data[i]; 
  }
  
  
  // Set the unweigthed cell averages and calculate the grand sum,
  // rows and columns averages, and factor interaction sum of squares
  
  T=0;
  Y=0;
  for (i=0;i<cls[0];i++) {
    for (j=0;j<cls[1];j++) {
      for (k=0;k<cls[2];k++) {
	T += cell_sum[i][j][k];
	Y += cell_ssq[i][j][k];
	col_sumA[j] += cell_sum[i][j][k];
	col_sumB[k] += cell_sum[i][j][k];
      }
    }
  }
  T=(T*T)/(cls[0]*cls[1]*cls[2]);
  
  // Set the sum of square for the 2 factor posibilities
  
  idx=0;
  for (i=0;i<cls[0];i++) {
    for (j=0;j<cls[1];j++) {
      for (k=0;k<cls[2];k++) {
	col_row_sumAS[idx] += cell_sum[i][j][k];
      }
      idx++;
    }
  }
  idx=0;
  for (j=0;j<cls[1];j++) {
    for (k=0;k<cls[2];k++) {
      for (i=0;i<cls[0];i++) {
	col_row_sumAB[idx] += cell_sum[i][j][k];
      }
      idx++;
    }
  }
  
  // Calculate the factors sum of squares
  
  A=0;
  for (i=0;i<cls[1];i++) {
    A += col_sumA[i]*col_sumA[i];
  }    
  A /= cls[0]*cls[2];
  
  B=0;
  for (i=0;i<cls[2];i++) {
    B += col_sumB[i]*col_sumB[i];
  }    
  B /= cls[0]*cls[1];
  
  AS=0;
  for (i=0;i<(cls[0]*cls[1]);i++) {
    AS += col_row_sumAS[i]*col_row_sumAS[i];
  }    
  AS /= cls[2];
  
  AB=0;
  for (i=0;i<(cls[1]*cls[2]);i++) {
    AB += col_row_sumAB[i]*col_row_sumAB[i];
  }
  AB /= cls[0];
  
  // Sum of Squares
  
  SSA = A-T+TINY;
  SSB = B-T+TINY;
  SSAB = AB-A-B+T+TINY;
  SSAS = AS-A+TINY;
  SSY = Y-AB-AS+A+TINY;
  SST = Y-T+TINY;
  
  // Degrees of freedom
  
  dfA = cls[1]-1+TINY;
  dfB = cls[2]-1+TINY;
  dfAB = dfA*dfB;
  dfAS = (cls[1]*(cls[0]-1))+TINY;
  dfY = dfAS*dfB;
  dfT = n-1;
  
  // Mean sum of Squares
  
  MSA=SSA/dfA;
  MSB=SSB/dfB;
  MSAB=SSAB/dfAB;
  MSAS=SSAS/dfAS;
  MSY=SSY/dfY;
  
  // F-ratios
  
  fA=MSA/MSAS;
  fB=MSB/MSY;
  fAB=MSAB/MSY;
  
  // Probabilities
  
  *probA=betai(0.5*dfAS,0.5*dfA,dfAS/(dfAS+dfA*(fA)));
  *probB=betai(0.5*dfY,0.5*dfB,dfY/(dfY+dfB*(fB)));
  *probAB=betai(0.5*dfY,0.5*dfAB,dfY/(dfY+dfAB*(fAB)));
  
  
  if (verbose || debug) {
    printf ("DATA\n");
    print_anova_data(data, n, factor, 3);
    printf ("SOURCE                               SS     df           MS         F    Pr(>F)\n");
    strcpy(str, *nameA);
    strcat(str, " (A)");
    printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSA, dfA, MSA, fA, *probA);
    strcpy(str, *nameB);
    strcat(str, " (B)");
    printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSB, dfB, MSB, fB, *probB);
    strcpy(str, "A*B");
    printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSAB, dfAB, MSAB, fAB, *probAB);
    strcpy(str, *nameS);
    strcat(str, " (S)/A");
    printf ("%-24s %14.2f %6i %12.2f\n", str, SSAS, dfAS, MSAS);
    strcpy(str, "B*S/A");
    printf ("%-24s %14.2f %6i %12.2f\n", str, SSY, dfY, MSY);
    strcpy(str, "Total");
    printf ("%-24s %14.2f %6i\n", str, SST, dfT);
  }
  
}

void
anova_two_way_within(double data[], int **factor, int n, int cls[],
		     double *probA, double *probB, double *probAB,
		     char **nameS, char **nameA, char **nameB)
{
  
  // This code can only support balanced within-subjects designs

  int i, j, k, idx;
  double T, A, B, S, AB, AS, BS, Y;
  double SST, SSA, SSB, SSS, SSAB, SSAS, SSBS, SSABS;
  double MSA, MSB, MSS, MSAB, MSAS, MSBS, MSABS;
  int dfT, dfA, dfB, dfS, dfAB, dfAS, dfBS, dfABS;
  double fA, fB, fAB;
  double col_sumS[cls[0]];
  double col_sumA[cls[1]];
  double col_sumB[cls[2]];
  double col_row_sumAS[cls[0]*cls[1]];
  double col_row_sumBS[cls[0]*cls[2]];
  double col_row_sumAB[cls[1]*cls[2]];
  double cell_sum[cls[0]][cls[1]][cls[2]];
  double cell_ssq[cls[0]][cls[1]][cls[2]];
  char	str[25];
  
  // Initialize the arrays otherwise weird things may happen
  
  for (i=0;i<cls[0];i++) {
    col_sumS[i] = 0;
  }
  for (i=0;i<cls[1];i++) {
    col_sumA[i] = 0;    
  }
  for (i=0;i<cls[2];i++) {
    col_sumB[i] = 0;    
  }
  for (i=0;i<(cls[0]*cls[1]);i++) {
    col_row_sumAS[i] = 0;
  }
  for (i=0;i<(cls[0]*cls[2]);i++) {
    col_row_sumBS[i] = 0;
  }
  for (i=0;i<(cls[1]*cls[2]);i++) {
    col_row_sumAB[i] = 0;
  }
  for (i=0;i<cls[0];i++) {
    for (j=0;j<cls[1];j++) {
      for (k=0;k<cls[2];k++) {
	cell_sum[i][j][k] = 0;
	cell_ssq[i][j][k] = 0;
      }
    }
  }
  
  // Calculate cell sums and cell sum of squares
  
  for (i=0;i<n;i++) {
    cell_ssq[factor[0][i]][factor[1][i]][factor[2][i]] += data[i]*data[i]; 
    cell_sum[factor[0][i]][factor[1][i]][factor[2][i]] += data[i]; 
  }
  
  
  // Set the unweigthed cell averages and calculate the grand sum,
  // rows and columns averages, and factor interaction sum of squares
  
  T=0;
  Y=0;
  for (i=0;i<cls[0];i++) {
    for (j=0;j<cls[1];j++) {
      for (k=0;k<cls[2];k++) {
	T += cell_sum[i][j][k];
	Y += cell_ssq[i][j][k];
	col_sumS[i] += cell_sum[i][j][k];
	col_sumA[j] += cell_sum[i][j][k];
	col_sumB[k] += cell_sum[i][j][k];
      }
    }
  }
  T=(T*T)/(cls[0]*cls[1]*cls[2]);
  
  // Set the sum of square for the 2 factor posibilities
  
  idx=0;
  for (i=0;i<cls[0];i++) {
    for (j=0;j<cls[1];j++) {
      for (k=0;k<cls[2];k++) {
	col_row_sumAS[idx] += cell_sum[i][j][k];
      }
      idx++;
    }
  }
  idx=0;
  for (i=0;i<cls[0];i++) {
    for (k=0;k<cls[2];k++) {
      for (j=0;j<cls[1];j++) {
	col_row_sumBS[idx] += cell_sum[i][j][k];
      }
      idx++;
    }
  }
  idx=0;
  for (j=0;j<cls[1];j++) {
    for (k=0;k<cls[2];k++) {
      for (i=0;i<cls[0];i++) {
	col_row_sumAB[idx] += cell_sum[i][j][k];
      }
      idx++;
    }
  }
  
  // Calculate the factors sum of squares
  
  S=0;
  for (i=0;i<cls[0];i++) {
    S += col_sumS[i]*col_sumS[i];
  }
  S /= cls[1]*cls[2];
  
  A=0;
  for (i=0;i<cls[1];i++) {
    A += col_sumA[i]*col_sumA[i];
  }    
  A /= cls[0]*cls[2];
  
  B=0;
  for (i=0;i<cls[2];i++) {
    B += col_sumB[i]*col_sumB[i];
  }    
  B /= cls[0]*cls[1];
  
  AS=0;
  for (i=0;i<(cls[0]*cls[1]);i++) {
    AS += col_row_sumAS[i]*col_row_sumAS[i];
  }    
  AS /= cls[2];
  
  BS=0;
  for (i=0;i<(cls[0]*cls[2]);i++) {
    BS += col_row_sumBS[i]*col_row_sumBS[i];
  }    
  BS /= cls[1];
  
  AB=0;
  for (i=0;i<(cls[1]*cls[2]);i++) {
    AB += col_row_sumAB[i]*col_row_sumAB[i];
  }
  AB /= cls[0];
  
  // Sum of Squares
  
  SSS = S-T+TINY;
  SSA = A-T+TINY;
  SSB = B-T+TINY;
  SSAS = AS-A-S+T+TINY;
  SSBS = BS-B-S+T+TINY;
  SSAB = AB-A-B+T+TINY;
  SSABS = Y-AS-BS-AB+A+B+S-T+TINY;
  SST = SSS+SSA+SSB+SSAS+SSBS+SSAB+SSABS+TINY;
  
  // Degrees of freedom
  
  dfS = cls[0]-1+TINY;
  dfA = cls[1]-1+TINY;
  dfB = cls[2]-1+TINY;
  dfAS = dfS*dfA;
  dfBS = dfS*dfB;
  dfAB = dfA*dfB;
  dfABS = dfS*dfA*dfB;
  dfT = n-1;
  
  // Mean sum of Squares
  
  MSS=SSS/dfS;
  MSA=SSA/dfA;
  MSB=SSB/dfB;
  MSAS=SSAS/dfAS;
  MSBS=SSBS/dfBS;
  MSAB=SSAB/dfAB;
  MSABS=SSABS/dfABS;
  
  // F-ratios
  
  fA=MSA/MSAS;
  fB=MSB/MSBS;
  fAB=MSAB/MSABS;

  // Probabilities
  
  *probA=betai(0.5*dfABS,0.5*dfA,dfABS/(dfABS+dfA*(fA)));
  *probB=betai(0.5*dfABS,0.5*dfB,dfABS/(dfABS+dfB*(fB)));
  *probAB=betai(0.5*dfABS,0.5*dfAB,dfABS/(dfABS+dfAB*(fAB)));
  
  
  if (verbose || debug) {
    printf ("DATA\n");
    print_anova_data(data, n, factor, 3);
    printf ("\nSOURCE                               SS     df           MS         F    Pr(>F)\n");
    strcpy(str, *nameA);
    strcat(str, " (A)");
    printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSA, dfA, MSA, fA, *probA);
    strcpy(str, *nameB);
    strcat(str, " (B)");
    printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSB, dfB, MSB, fB, *probB);
    strcpy(str, "A*B");
    printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSAB, dfAB, MSAB, fAB, *probAB);
    strcpy(str, *nameS);
    strcat(str, " (S)");
    printf ("%-24s %14.2f %6i %12.2f\n", str, SSS, dfS, MSS);
    strcpy(str, "A*S");
    printf ("%-24s %14.2f %6i %12.2f\n", str, SSAS, dfAS, MSAS);
    strcpy(str, "B*S");
    printf ("%-24s %14.2f %6i %12.2f\n", str, SSBS, dfBS, MSBS);
    strcpy(str, "A*B*S");
    printf ("%-24s %14.2f %6i %12.2f\n", str, SSABS, dfABS, MSABS);
    strcpy(str, "Total");
    printf ("%-24s %14.2f %6i\n", str, SST, dfT);
  }
  
}

void
anova_three_way(double data[], int **factor, int n_tot, int cls[], int inter, 
		double *probA, double *probB, double *probC, 
		double *probAB, double *probAC, double *probBC, double *probABC,
		char **nameA, char **nameB, char **nameC)
{
  
  // This code utilizes the traditional unweighted means

  int i, j, k, idx;
  double T, A, B, C, AB, AC, BC, ABC, Y, n;
  double SST, SSA, SSB, SSC, SSAB, SSAC, SSBC, SSABC, SSY;
  double MSA, MSB, MSC, MSAB, MSAC, MSBC, MSABC, MSY;
  int dfT, dfA, dfB, dfC, dfAB, dfAC, dfBC, dfABC, dfY;
  double fA, fB, fC, fAB, fAC, fBC, fABC;
  double col_sumA[cls[0]];
  double col_sumB[cls[1]];
  double col_sumC[cls[2]];
  double col_row_sumAB[cls[0]*cls[1]];
  double col_row_sumAC[cls[0]*cls[2]];
  double col_row_sumBC[cls[1]*cls[2]];
  double cell_sum[cls[0]][cls[1]][cls[2]];
  double cell_ssq[cls[0]][cls[1]][cls[2]];
  double cell_avg[cls[0]][cls[1]][cls[2]];
  double cell_uav[cls[0]][cls[1]][cls[2]];
  double cell_mem[cls[0]][cls[1]][cls[2]];
  char	str[25];
  
  // Initialize the arrays otherwise weird things may happen
  
  for (i=0;i<cls[0];i++) {
    col_sumA[i] = 0;
  }
  for (i=0;i<cls[1];i++) {
    col_sumB[i] = 0;    
  }
  for (i=0;i<cls[2];i++) {
    col_sumC[i] = 0;    
  }
  for (i=0;i<(cls[0]*cls[1]);i++) {
    col_row_sumAB[i] = 0;
  }
  for (i=0;i<(cls[0]*cls[2]);i++) {
    col_row_sumAC[i] = 0;
  }
  for (i=0;i<(cls[1]*cls[2]);i++) {
    col_row_sumBC[i] = 0;
  }
  for (i=0;i<cls[0];i++) {
    for (j=0;j<cls[1];j++) {
      for (k=0;k<cls[2];k++) {
	cell_sum[i][j][k] = 0;
	cell_ssq[i][j][k] = 0;
	cell_avg[i][j][k] = 0;
	cell_uav[i][j][k] = 0;
	cell_mem[i][j][k] = 0;
      }
    }
  }
  
  // Calculate cell sums and cell sum of squares

  for (i=0;i<n_tot;i++) {
    cell_ssq[factor[0][i]][factor[1][i]][factor[2][i]] += data[i]*data[i]; 
    cell_sum[factor[0][i]][factor[1][i]][factor[2][i]] += data[i]; 
    cell_mem[factor[0][i]][factor[1][i]][factor[2][i]]++;
  }
  
  // Calculate cell averages, harmonic mean, within groups degrees of freedom
  // and sum of squares for unbalanced designs
  
  n=0;
  Y=0;
  dfY=0;
  for (i=0;i<cls[0];i++) {
    for (j=0;j<cls[1];j++) {
      for (k=0;k<cls[2];k++) {
	cell_avg[i][j][k]=cell_sum[i][j][k]/cell_mem[i][j][k];
	n += (1/cell_mem[i][j][k]);
	Y += cell_ssq[i][j][k]-(cell_sum[i][j][k]*cell_sum[i][j][k]/cell_mem[i][j][k]);
	dfY += (cell_mem[i][j][k]-1);
      }
    }
  }
  n=cls[0]*cls[1]*cls[2]/n;
  
  // Set the unweigthed cell averages and calculate the grand sum,
  // rows and columns averages, and factor interaction sum of squares
  
  T=0;
  ABC=0;
  for (i=0;i<cls[0];i++) {
    for (j=0;j<cls[1];j++) {
      for (k=0;k<cls[2];k++) {
	cell_uav[i][j][k]=cell_avg[i][j][k]*n;
	T += cell_uav[i][j][k];
	ABC += cell_uav[i][j][k]*cell_uav[i][j][k];
	col_sumA[i] += cell_uav[i][j][k];
	col_sumB[j] += cell_uav[i][j][k];
	col_sumC[k] += cell_uav[i][j][k];
      }
    }
  }
  T=(T*T)/(cls[0]*cls[1]*cls[2]*n);
  ABC /= n;
  
  // Set the sum of square for the 2 factor posibilities
  
  idx=0;
  for (i=0;i<cls[0];i++) {
    for (j=0;j<cls[1];j++) {
      for (k=0;k<cls[2];k++) {
	col_row_sumAB[idx] += cell_uav[i][j][k];
      }
      idx++;
    }
  }
  idx=0;
  for (i=0;i<cls[0];i++) {
    for (k=0;k<cls[2];k++) {
      for (j=0;j<cls[1];j++) {
	col_row_sumAC[idx] += cell_uav[i][j][k];
      }
      idx++;
    }
  }
  idx=0;
  for (j=0;j<cls[1];j++) {
    for (k=0;k<cls[2];k++) {
      for (i=0;i<cls[0];i++) {
	col_row_sumBC[idx] += cell_uav[i][j][k];
      }
      idx++;
    }
  }
  
  // Calculate the factors sum of squares
  
  A=0;
  for (i=0;i<cls[0];i++) {
    A += col_sumA[i]*col_sumA[i];
  }
  A /= cls[1]*cls[2]*n;
  
  B=0;
  for (i=0;i<cls[1];i++) {
    B += col_sumB[i]*col_sumB[i];
  }    
  B /= cls[0]*cls[2]*n;
  
  C=0;
  for (i=0;i<cls[2];i++) {
    C += col_sumC[i]*col_sumC[i];
  }    
  C /= cls[0]*cls[1]*n;
  
  AB=0;
  for (i=0;i<(cls[0]*cls[1]);i++) {
    AB += col_row_sumAB[i]*col_row_sumAB[i];
  }
  AB /= cls[2]*n;
  
  AC=0;
  for (i=0;i<(cls[0]*cls[2]);i++) {
    AC += col_row_sumAC[i]*col_row_sumAC[i];
  }    
  AC /= cls[1]*n;
  
  BC=0;
  for (i=0;i<(cls[1]*cls[2]);i++) {
    BC += col_row_sumBC[i]*col_row_sumBC[i];
  }    
  BC /= cls[0]*n;
  
  // Sum of Squares
  
  SSA = A-T+TINY;
  SSB = B-T+TINY;
  SSC = C-T+TINY;
  SSAB = AB-A-B+T+TINY;
  SSAC = AC-A-C+T+TINY;
  SSBC = BC-B-C+T+TINY;
  SSABC = ABC-AB-AC-BC+A+B+C-T+TINY;
  if (! inter) {
    SSY = Y+SSAB+SSAC+SSBC+SSABC+TINY;
    SST = SSA+SSB+SSC+SSY+TINY;
  } else {
    SSY = Y+TINY;
    SST = SSA+SSB+SSC+SSAB+SSAC+SSBC+SSABC+SSY+TINY;
  }
  
  // Degrees of freedom
  
  dfA = cls[0]-1+TINY;
  dfB = cls[1]-1+TINY;
  dfC = cls[2]-1+TINY;
  dfAB = dfA*dfB;
  dfAC = dfA*dfC;
  dfBC = dfB*dfC;
  dfABC = dfA*dfB*dfC;
  dfT = n_tot-1;
  if (! inter) dfY=(n_tot-dfA-dfB-dfC-1)+TINY;
  
  // Mean sum of Squares
  
  MSA=SSA/dfA;
  MSB=SSB/dfB;
  MSC=SSC/dfC;
  MSAB=SSAB/dfAB;
  MSAC=SSAC/dfAC;
  MSBC=SSBC/dfBC;
  MSABC=SSABC/dfABC;
  MSY=SSY/dfY;
  
  // F-ratios
  
  fA=MSA/MSY;
  fB=MSB/MSY;
  fC=MSC/MSY;
  fAB=MSAB/MSY;
  fAC=MSAC/MSY;
  fBC=MSBC/MSY;
  fABC=MSABC/MSY;
  
  // Probabilities
  
  *probA=betai(0.5*dfY,0.5*dfA,dfY/(dfY+dfA*(fA)));
  *probB=betai(0.5*dfY,0.5*dfB,dfY/(dfY+dfB*(fB)));
  *probC=betai(0.5*dfY,0.5*dfC,dfY/(dfY+dfC*(fC)));
  *probAB=betai(0.5*dfY,0.5*dfAB,dfY/(dfY+dfAB*(fAB)));
  *probAC=betai(0.5*dfY,0.5*dfAC,dfY/(dfY+dfAC*(fAC)));
  *probBC=betai(0.5*dfY,0.5*dfBC,dfY/(dfY+dfBC*(fBC)));
  *probABC=betai(0.5*dfY,0.5*dfABC,dfY/(dfY+dfABC*(fABC)));
  
  if (verbose || debug) {
    printf ("DATA\n");
    print_anova_data(data, n_tot, factor, 3);
    printf ("\nSOURCE                               SS     df           MS         F    Pr(>F)\n");
    strcpy(str, *nameA);
    strcat(str, " (A)");
    printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSA, dfA, MSA, fA, *probA);
    strcpy(str, *nameB);
    strcat(str, " (B)");
    printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSB, dfB, MSB, fB, *probB);
    strcpy(str, *nameC);
    strcat(str, " (C)");
    printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSC, dfC, MSC, fC, *probC);
    if (inter) {
      strcpy(str, "A*B");
      printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSAB, dfAB, MSAB, fAB, *probAB);
      strcpy(str, "A*C");
      printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSAC, dfAC, MSAC, fAC, *probAC);
      strcpy(str, "B*C");
      printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSBC, dfBC, MSBC, fBC, *probBC);
      strcpy(str, "A*B*C");
      printf ("%-24s %14.2f %6i %12.2f %9.3f  %7.4g\n", str, SSABC, dfABC, MSABC, fABC, *probABC);
    } else {
      *probAB=1.0;
    }
    strcpy(str, "Within (error)");
    printf ("%-24s %14.2f %6i %12.2f\n", str, SSY, dfY, MSY);
    strcpy(str, "Total");
    printf ("%-24s %14.2f %6i\n", str, SST, dfT);
  }
  
}

int
main(int argc,
     char **argv)

{

  char *data_file, *factor_file, *output_file;
  int  samples, variables, inter, permuts, lines;
  enum TT tt;
  XDATA *xdata;
  FACTOR *factor;
  RESULT *res;
  
  verbose = 0;

  debug = 0;

  parse_arguments(argc, argv, &data_file, &factor_file, &output_file,
		  &samples, &variables, &inter, &tt, &permuts);
  
  if (samples == 0) {
    samples = count_fields(factor_file);
  }
  
  lines = count_lines(factor_file);
  lines--;

  factor = read_factor(factor_file, lines, samples);
  
  if (variables == 0) {
    variables = count_lines(data_file);
    variables--;
  }
  
  xdata = read_xdata(data_file, samples, variables);
  
  if (lines == 1) {
    res = run_one_way_anova(xdata, factor, permuts);
  } else {
    res = run_multi_factorial_anova(xdata, factor, tt, inter, permuts);
  }

  print_results(output_file, res, 1, res->a-1);   

  free(data_file);
  free(factor_file);
  free(output_file);
  free(xdata);
  free(factor);
  free(res);

  exit(0);

}
