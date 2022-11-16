/*********************************************************************
 Copyright 2002 Bristol-Myers Squibb. All rights reserved.
 
   Author: Isaac M. Neuhaus and Robert E. Bruccoleri
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Many thanks to Mark Friedrichs, Carlod Rios and John Hinsdale
   for their help during debugging this code.

 $Id: nna.c,v 1.21 2008/09/11 17:39:47 neuhausi Exp $
**********************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <malloc.h>
#define EXTERN
#include "nna.h"

void
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **profile_file,
		char **output_file,
		int *samples,
		int *variables,
		enum METRIC_ENUM *metrics,
		double *order,
		char ***obj_names,
		int *transpose,
		int *inverse,
		int *n_objects)
{
  
  int c, arg_len, max_objs;
  extern char *optarg;
  int errflg = 0;
  
  *data_file = NULL;
  
  *profile_file = NULL;
  
  *output_file = NULL;
  
  *samples = 0;
  
  *variables = 0;
  
  *metrics = EUCLIDEAN;
  
  *order = 2.0;
  
  /* Maximum of 100 objects to make a profile */
  max_objs = 100;
  
  *obj_names = ARRAY_ALLOC(max_objs, char *);
  
  *transpose = 0;
  
  *inverse = 0;
  
  *n_objects = 0;
  
  while ((c = getopt(argc, argv, "d:p:o:s:g:m:r:n:tivD")) != EOF)
    switch (c) {
    case 'd':
      *data_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*data_file, optarg);
      break;
    case 'p':
      *profile_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*profile_file, optarg);
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
    case 'm':
      if (strcmp(optarg, "euclidean") == 0) {
	*metrics = EUCLIDEAN;
      } else if (strcmp(optarg, "minkowski") == 0) {
	*metrics = MINKOWSKI;
      } else if (strcmp(optarg, "manhattan") == 0) {
	*metrics = MANHATTAN;
      } else if (strcmp(optarg, "maximum") == 0) {
	*metrics = MAXIMUM;
      } else if (strcmp(optarg, "pearson") == 0) {
	*metrics = PEARSON;
      } else if (strcmp(optarg, "spearman") == 0) {
	*metrics = SPEARMAN;
      } else if (strcmp(optarg, "lspearman") == 0) {
	*metrics = LSPEARMAN;
      } else {
	errflg++;
      }
      break;
    case 'r':
      *order = atof(optarg);
      break;
    case 'n':
      if (*n_objects <= max_objs) {
	arg_len = strlen(optarg);
	(*obj_names)[*n_objects] = ARRAY_ALLOC((arg_len+1), char);
	strcpy((*obj_names)[*n_objects], optarg);
	*n_objects += 1;
      }
      break;
    case 't':
      *transpose = 1;
      break;
    case 'i':
      *inverse = 1;
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
    }
  
  if (*data_file == NULL ) { 
    errflg++;
  }
  
  if (*n_objects == 0 && *profile_file == NULL) { 
    errflg++;
  }
  
  if (errflg) {
    printf("\nUsage: nna <options> (options in square brackets are optional)\n\n");
    printf("The nearest neighbor algorithm, nna, can be used to find similar patterns in\n");
    printf("the rows of a specified data file. The pattern can be obtained from either\n");
    printf("any row in the data file or defined in a different profile file.\n");
    printf(" -d  <string>    data file. A tab delimeted file with the following format.\n");
    printf("                 A header line with sample names. Lines with expression values\n");
    printf("                 preceded by the variable name at the begining of the line.\n");
    printf(" -p  <string>    The file containing the profile pattern\n");
    printf(" -o  [<string>]  The output file (The output file name will be the string given\n");
    printf("                 whit the extension '.nna'. If no file is specified it defaults\n");
    printf("                 to stdout.\n");
    printf(" -s  [<integer>] number of samples. It should correspond to the number of names\n");
    printf("                 in the first line of the data file.\n");
    printf(" -g  [<integer>] number of variables. It should correspond to the number of lines + 1\n");
    printf("                 in the data file.\n");
    printf(" -m  [string]    The metrics to use for the nna. The available metrics are\n");
    printf("                 euclidean, minkowski, manhattan, maximum, pearson, spearman or\n");
    printf("                 lspearman (log spearman).\n");
    printf("                 The default is euclidean\n");
    printf(" -r  [float]     If the minkowski metric is selected, a power may be specified.\n");
    printf("                 The default is 2.0 (which is equivalent to the euclidean distance)\n");
    printf(" -n  <string>    The profile name\n");
    printf(" -t  [switch]    Transpose data file. The default is not to transpose\n");
    printf(" -i  [switch]    Find similar or dissimilar objects. The default\n");
    printf("                 is off to find similar objects\n\n");
    printf("Examples\n");
    printf("\n\tnna -d ./t/nna.dat -n M55593_at\n\n");
    exit(0);
  }
    
}

VAR_VECTOR *
create_var_pattern(int n_objects, int *obj_numbers, XDATA *xdata)
{

  int i, j, n;
  int vector_size = xdata->s;
  double mean, stdev, *val;
  VECTOR *v;
  VAR_VECTOR *vv;
  
  vv = TYPE_ALLOC(VAR_VECTOR);
  vv->mean = ARRAY_ALLOC(vector_size, double);
  vv->stdev = ARRAY_ALLOC(vector_size, double);
  vv->n = vector_size;
  
  for (i = 0; i < vector_size; i++) {
    mean = 0.0;
    stdev = 0.0;
    n = 0;
    for (j = 0; j < n_objects; j++) {
      v = xdata->vectors[obj_numbers[j]];
      val = v->data;
      if (! isnan(val[i])) {
	mean += val[i];
        stdev += val[i] * val[i];
	n++;
      }
    }
    if (n == 0) {
      vv->mean[i] = 0.0;
      vv->stdev[i] = 0.0;
    } else {
      mean /= n;
      vv->mean[i] = mean;
      vv->stdev[i] = sqrt(MAX(0.0, stdev / n - mean * mean));
    }
  }

  return vv;

}

int
search_obj(XDATA *xdata,
	   char *obj_name)
{
  int i;
  
  for (i = 0; i < xdata->n; i++) {
    if (strcmp(obj_name, xdata->variables[i]) == 0) {
      return i;
    }
  }
  
  return (int) -1;

}
	
void
calculate_nna(XDATA *xdata,
	      VECTOR *pattern,
	      char *name,
	      char *output_file,
	      enum METRIC_ENUM metrics,
	      double power,
	      int inverse)
{

  int i;
  double nna[xdata->n], nn;   
  double sum = 0.0;
  double sum2 = 0.0;
  
  if (inverse == 1 ) {
    for (i = 0; i < xdata->n; i++) {
      if (metrics == EUCLIDEAN) {
	nn = ieuclidean(pattern,
			xdata->vectors[i]);
      } else if (metrics == PEARSON) {
	nn = ipearson(pattern,
		      xdata->vectors[i]);
      } else if (metrics == SPEARMAN) {
	nn = ispearman(pattern,
		       xdata->vectors[i]);
      } else if (metrics == LSPEARMAN) {
	nn = ilspearman(pattern,
		        xdata->vectors[i]);
      } else if (metrics == MAXIMUM) {
	nn = imaximum(pattern,
		      xdata->vectors[i]);
      } else if (metrics == MANHATTAN) {
	nn = imanhattan(pattern,
			xdata->vectors[i]);
      } else {
	nn = iminkowski(pattern,
			xdata->vectors[i],
			power);
      }
      nna[i] = nn;
      sum += nn;
      sum2 += nn * nn;
    }
  } else {
    for (i = 0; i < xdata->n; i++) {
      if (metrics == EUCLIDEAN) {
	nn = euclidean(pattern,
		       xdata->vectors[i]);
      } else if (metrics == PEARSON) {
	nn = pearson(pattern,
		     xdata->vectors[i]);
      } else if (metrics == SPEARMAN) {
	nn = spearman(pattern,
		      xdata->vectors[i]);
      } else if (metrics == LSPEARMAN) {
	nn = lspearman(pattern,
		       xdata->vectors[i]);
      } else if (metrics == MAXIMUM) {
	nn = maximum(pattern,
		     xdata->vectors[i]);
      } else if (metrics == MANHATTAN) {
	nn = manhattan(pattern,
		       xdata->vectors[i]);
      } else {
	nn = minkowski(pattern,
		       xdata->vectors[i],
		       power);
      }
      nna[i] = nn;
      sum += nn;
      sum2 += nn * nn;
    }    
  }

  print_nna_results(xdata, name, output_file, metrics, sum, sum2, nna);

}

void
calculate_nna_var(XDATA *xdata,
		  VAR_VECTOR *var_pattern,
		  char *name,
		  char *output_file,
		  enum METRIC_ENUM metrics,
		  double power,
		  int inverse)
{

  int i;
  double nna[xdata->n], nn;   
  double sum = 0.0;
  double sum2 = 0.0;
  
  if (inverse == 1 ) {
    for (i = 0; i < xdata->n; i++) {
      if (metrics == EUCLIDEAN) {
	nn = var_ieuclidean(var_pattern,
			    xdata->vectors[i]);
      } else if (metrics == PEARSON) {
	nn = var_ipearson(var_pattern,
			  xdata->vectors[i]);
      } else if (metrics == SPEARMAN) {
	nn = var_ispearman(var_pattern,
			   xdata->vectors[i]);
      } else if (metrics == LSPEARMAN) {
	nn = var_ilspearman(var_pattern,
			    xdata->vectors[i]);
      } else if (metrics == MAXIMUM) {
	nn = var_imaximum(var_pattern,
			  xdata->vectors[i]);
      } else if (metrics == MANHATTAN) {
	nn = var_imanhattan(var_pattern,
			    xdata->vectors[i]);
      } else {
	nn = var_iminkowski(var_pattern,
			    xdata->vectors[i],
			    power);
      }
      nna[i] = nn;
      sum += nn;
      sum2 += nn * nn;
    }
  } else {
    for (i = 0; i < xdata->n; i++) {
      if (metrics == EUCLIDEAN) {
	nn = var_euclidean(var_pattern,
			   xdata->vectors[i]);
      } else if (metrics == PEARSON) {
	nn = var_pearson(var_pattern,
			 xdata->vectors[i]);
      } else if (metrics == SPEARMAN) {
	nn = var_spearman(var_pattern,
			  xdata->vectors[i]);
      } else if (metrics == LSPEARMAN) {
	nn = var_lspearman(var_pattern,
			   xdata->vectors[i]);
      } else if (metrics == MAXIMUM) {
	nn = var_maximum(var_pattern,
			 xdata->vectors[i]);
      } else if (metrics == MANHATTAN) {
	nn = var_manhattan(var_pattern,
			   xdata->vectors[i]);
      } else {
	nn = var_minkowski(var_pattern,
			   xdata->vectors[i],
			   power);
      }
      nna[i] = nn;
      sum += nn;
      sum2 += nn * nn;
    }    
  }

  print_nna_results(xdata, name, output_file, metrics, sum, sum2, nna);

}

void
print_nna_results(XDATA *xdata,
		  char *name,
		  char *output_file,
		  enum METRIC_ENUM metrics,
		  double sum,
		  double sum2,
		  double *nna)
{

  int i, j, perm[xdata->n];
  char *filename;
  FILE *out;
  double mean, stdv;

  mean = sum / xdata->n;
  stdv = sqrt(MAX(0.0, sum2 / xdata->n - mean * mean));
  sort_by_index(xdata->n, nna, perm);

  if (output_file == NULL) {

    printf("NN for %s\t%5.4f\t%5.4f\t%d\t%d\n",
	   name,
	   mean,
	   stdv,
	   xdata->n,
	   xdata->s);

    if (metrics == EUCLIDEAN || metrics == MINKOWSKI ||
	metrics == MANHATTAN || metrics == MAXIMUM) {
      for (i = 0; i < xdata->n; i++) {
	j = perm[i];
	printf("%s\t%5.4f\t%5.4f\n",
	       xdata->variables[j],
	       nna[j],
	       stdv > 0.0 ? (nna[j] - mean) / stdv : 0.0);
      }
    } else {
      for (i = xdata->n - 1; i >= 0; i--) {
	j = perm[i];
	printf("%s\t%5.8f\t%5.8f\n",
	       xdata->variables[j],
	       nna[j],
	       stdv > 0.0 ? (nna[j] - mean) / stdv : 0.0);
      }
    }

  } else {

    filename = ARRAY_ALLOC((strlen(output_file) + strlen(".nna") + 1), char);
    strcpy(filename, output_file);
    strcat(filename, ".nna");
  
    if ((out = fopen(filename, "w")) == NULL) {
      exit_on_error("Error writing output file");
    }

    fprintf(out,
	    "NN for %s\t%5.4f\t%5.4f\t%d\t%d\n",
	    name,
	    mean,
	    stdv,
	    xdata->n,
	    xdata->s);

    if (metrics == EUCLIDEAN || metrics == MINKOWSKI ||
	metrics == MANHATTAN || metrics == MAXIMUM) {
      for (i = 0; i < xdata->n; i++) {
	j = perm[i];
	fprintf(out,
		"%s\t%5.4f\t%5.4f\n",
		xdata->variables[j],
		nna[j],
		stdv > 0.0 ? (nna[j] - mean) / stdv : 0.0);
      }
    } else {
      for (i = xdata->n - 1; i >= 0; i--) {
	j = perm[i];
	fprintf(out,
		"%s\t%5.8f\t%5.8f\n",
		xdata->variables[j],
		nna[j],
		stdv > 0.0 ? (nna[j] - mean) / stdv : 0.0);
      }
    }

    fclose(out);
    free(filename);

  }

}

int
main(int argc,
     char **argv)

{
  
  char *data_file, *profile_file, *output_file;
  int samples, variables, binary;
  enum METRIC_ENUM metrics;
  double order;
  char *name, **obj_names, tmp_name[10000]; 
  int  transpose, inverse, n_objects;
  int i, *obj_numbers;
  int obj_number, len;
  
  XDATA *xdata;
  VECTOR *pattern;
  VAR_VECTOR *var_pattern;
  
  verbose = 0;
  debug = 0;
  binary = 0;
  
  parse_arguments(argc, argv, &data_file, &profile_file, &output_file, 
		  &samples, &variables, &metrics, &order, &obj_names, 
		  &transpose, &inverse, &n_objects);
  
  binary = is_binary_file(data_file);
  
  if (binary > 0) {
    xdata = read_binary_xdata(data_file);
    variables = xdata->n;
    samples= xdata->s; 
  } else {
    if (samples == 0) {
      samples = count_fields(data_file);
    }
    if (variables == 0) {
      variables = count_lines(data_file);
      variables--;
    }
    xdata = read_xdata(data_file, samples, variables);
  }
  
  if (transpose == 1) {
    xdata = transpose_xdata(xdata);
  }
  
  if (n_objects == 0 && profile_file != NULL){

    pattern = read_pattern(profile_file, samples);
    len = strlen("Pattern");
    name = ARRAY_ALLOC(len, char);
    strcpy(name, "Patern");
    name[len] = '\0';
    calculate_nna(xdata, pattern, name, output_file, metrics, order, inverse);
    free(profile_file);
    free(pattern);

  } else if (n_objects == 1) {

    obj_number = search_obj(xdata, obj_names[0]);
    if (obj_number < 0) {
      fprintf(stderr, "Unable to find %s in data file %s\n", obj_names[0], data_file);
      exit(1);
    }
    pattern = xdata->vectors[obj_number];
    len = strlen(xdata->variables[obj_number]);
    name = ARRAY_ALLOC(len, char);
    strcpy(name, xdata->variables[obj_number]);
    calculate_nna(xdata, pattern, name, output_file, metrics, order, inverse);
    free(pattern);
    free(obj_names);

  } else if (n_objects > 1){

    len = 0;
    strcpy(tmp_name, "");
    obj_numbers = ARRAY_ALLOC(n_objects, int);
    for (i = 0; i < n_objects; i++) {
      obj_number = search_obj(xdata, obj_names[i]);
      strcat(tmp_name, xdata->variables[obj_number]); 
      len += strlen(xdata->variables[obj_number]);
      if (i < n_objects -1) {
	strcat(tmp_name, "+"); 
	len += 1;
      }
      obj_numbers[i] = obj_number;
      if (obj_number < 0) {
	fprintf(stderr,	"Unable to find %s in data file %s\n", obj_names[i], data_file);
	exit(1);
      }
    }
    var_pattern = create_var_pattern(n_objects, obj_numbers, xdata);
    name = ARRAY_ALLOC(len, char);
    strcpy(name, tmp_name);
    calculate_nna_var(xdata, var_pattern, name, output_file, metrics, order, inverse);
    free(var_pattern);
    free(obj_names);
    free(obj_numbers);

  }

  free(data_file);
  free(name);
  free(xdata);

  if (output_file != NULL) {
    free(output_file);
  }

  exit(0);
 
}
