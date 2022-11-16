/*********************************************************************
 Copyright 2006 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: ks.c,v 1.11 2008/03/24 13:55:28 neuhausi Exp $
**********************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define EXTERN
#include "ks.h"

void
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **output_file,
		int *wells,
		int *features,
		char ***reference,
		int *n_reference,
		char **feature_name)
{
  
  int c, arg_len, max;
  extern char *optarg;
  int errflg = 0;
  
  *data_file = NULL;
  
  *output_file = NULL;

  *feature_name = NULL;
  
  *wells = 0;
  
  *features = 0;
  
  *n_reference = 0;

  max = 100;

  *reference = ARRAY_ALLOC(max, char *);
  
  while ((c = getopt(argc, argv, "d:o:w:f:r:n:svD")) != EOF)
    switch (c) {
    case 'd':
      *data_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*data_file, optarg);
      break;
    case 'o':
      *output_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*output_file, optarg);
      break;
    case 'w':
      *wells = atoi(optarg);
      break;
    case 'f':
      *features = atoi(optarg);
      break;
    case 'r':
      if (*n_reference <= max) {
	arg_len = strlen(optarg);
	(*reference)[*n_reference] = ARRAY_ALLOC((arg_len+1), char);
	strcpy((*reference)[*n_reference], optarg);
	*n_reference += 1;
      }
      break;
    case 'n':
      *feature_name = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*feature_name, optarg);
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

  if (*data_file == NULL || *n_reference == 0) { 
    errflg++;
  }

  if (errflg) {
    printf("\nKS: Kolmogorov-Smirnov algorithm for Cellomics data\n");
    printf("\nUsage: ks <options> (options in square brackets are optional)\n\n");
    printf(" -d  <string>    data file. A tab delimeted file with the following format.\n");
    printf("                 A header line with feature names. Lines with intensity values\n");
    printf("                 preceded by the well identifier at the begining of the line.\n");
    printf("                 This identifier should be in the following format A1.xxx.xxx,\n");
    printf("                 where A1 is the well Id and the rest in discarded info \n");
    printf(" -o  [<string>]  output file. If no file is specified it defaults to stdout.\n");
    printf(" -w  [<integer>] number of wells. It should correspond to the number of wells\n");
    printf("                 the data will be summarized.\n");
    printf(" -f  [<integer>] number of features. It should correspond to the number of fields\n");
    printf("                 in the header in the data file.\n");
    printf(" -r  <string>    well reference (must match one or more names of the wells in the\n");
    printf("                 data file. Multiple reference wells are entered with multiple -r\n");
    printf("                 in the command\n");
    printf(" -n  <string>    feature name. Run ks only for the feature name\n");
    printf("Examples\n");
    printf("\tks -d ./t/ks.dat -r A0\n");
    printf("\n");
    exit(0);
  }
   
}

void
kstwo(double data1[], int n1, double data2[], int n2, double *d)
{

  int j1=0,j2=0;
  double d1,d2,dt,en1,en2,fn1=0.0,fn2=0.0,diff1=0.0,diff2=0.0;

  simple_sort(n2,data2);
  en1=n1;
  en2=n2;
  *d=0.0;

  // Check that not all the data has the same values

  diff1 = data1[n1-1] - data1[0];
  diff2 = data2[n2-1] - data2[0];


  if ((diff1 == 0.0) && (diff2 == 0.0)) {
    if (fabs(diff1 - diff2) > 0.0) {
      *d=1.0;
    }
    return;
  }

  while (j1 <= n1 && j2 <= n2) {
    d1 = data1[j1];
    d2 = data2[j2];
    if (d1 <= d2) {
      fn1 = j1++/en1;
    }
    if (d2 <= d1) {
      fn2 = j2++/en2;
    }
    dt = fabs(fn2 - fn1);
    if (dt > fabs(*d)) {
      *d = fn1 - fn2;
    }
  }

}

RESULT *
set_result(KSDATA *ksdata,
	   char *feature_name)
{

  int i, len;

  RESULT *res;
  VECTOR *vr;

  res = TYPE_ALLOC(RESULT);
  if (feature_name == NULL) {
    res->n = ksdata->f;
  } else {
    res->n = 1;
  }
  res->a = ksdata->w;
  res->attributes = ARRAY_ALLOC(res->a, char *);
  res->variables = ARRAY_ALLOC(res->n, char *);
  res->vectors = ARRAY_ALLOC(res->n, VECTOR *);
  for (i=0;i<res->n;i++) {
    vr = res->vectors[i] = TYPE_ALLOC(VECTOR);
    vr->data = ARRAY_ALLOC(res->a, double);
    vr->n = res->a;
  }
  for (i=0;i<res->a;i++) {
    len = strlen(ksdata->wells[i]);
    res->attributes[i] = ARRAY_ALLOC((len+1), char);  
    strncpy(res->attributes[i], ksdata->wells[i], len);
    res->attributes[i][len] = '\0';
  }
  if (feature_name == NULL) {
    for (i=0;i<res->n;i++) {
      len = strlen(ksdata->features[i]);
      res->variables[i] = ARRAY_ALLOC((len+1), char);  
      strncpy(res->variables[i], ksdata->features[i], len);
      res->variables[i][len] = '\0';
    }
  } else {
    len = strlen(feature_name);
    res->variables[0] = ARRAY_ALLOC((len+1), char);  
    strncpy(res->variables[0], feature_name, len);
    res->variables[0][len] = '\0';
  }

  return res;

}

RESULT *
calculate_ks (KSDATA *ksdata,
	      char **reference,
	      int n_reference,
	      char *feature_name)
{

  int i, j, k, c, cc, idxs[n_reference];
  double *ref, test[100000], d;
  RESULT *res;
  VECTOR *vr;

  res = set_result(ksdata, feature_name);

  // Set the indeces for the reference wells
  // and allocate memory for them
  c = 0;
  for (i = 0; i < ksdata->w; i++) {
    for (j = 0; j < n_reference; j++) {
      if (!strcmp(ksdata->wells[i], reference[j])) {
	idxs[j] = i;
	c += ksdata->c[i];
	break;
      }
    }
  }

  ref  = ARRAY_ALLOC(c, double);

  // Process the data
  for (i = 0; i < ksdata->f; i++) {
    if (feature_name == NULL) {
      vr = res->vectors[i];
    } else {
      if (strcmp(feature_name, ksdata->features[i]) == 0) {
	vr = res->vectors[0];
      } else {
	continue;
      }
    }
    // Get the reference wells
    c = 0;
    for (j = 0; j < n_reference; j++) {
      for (k = 0; k < ksdata->c[idxs[j]]; k++) {
	if (! isnan(ksdata->data[i][idxs[j]][k])) {
	  ref[c] = ksdata->data[i][idxs[j]][k];
	  c++;
	}
      }
    }
    simple_sort(c, ref);
    for (j = 0; j < ksdata->w; j++) {
      cc = 0;
      for (k = 0; k < ksdata->c[j]; k++) {
	if (! isnan(ksdata->data[i][j][k])) {
	  test[cc] = ksdata->data[i][j][k];
	  cc++;
	}
      }
      kstwo(ref, c, test, cc, &d);
      vr->data[j] = d;
    }
  }

  return res;

}

int 
main(int argc,
     char **argv)

{

  char *data_file, *output_file;
  int  wells, features, n_reference;
  char **reference, *feature_name;
  KSDATA *ksdata;
  RESULT *res;
  
  verbose = 0;

  debug = 0;

  parse_arguments(argc, argv,
		  &data_file,
		  &output_file,
		  &wells, &features,
		  &reference, &n_reference,
		  &feature_name);
    
  if (wells == 0) {
    wells = count_wells_ks(data_file);
  }

  if (features == 0) {
    features = count_fields(data_file);
  }

  ksdata = read_ksdata(data_file, features, wells);

  res = calculate_ks(ksdata, reference, n_reference, feature_name);

  print_results(output_file, res, 0 , 0);
  
  free(data_file);
  free(output_file);
  free(feature_name);
  free(reference);
  free(ksdata);
  free(res);

  exit(0);

}
