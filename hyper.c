/*********************************************************************
 Copyright 2007 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus and Roumyana Yordanova
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net
           roumyana.yordanova@bms.com or roumyana.yordanova@gmail.com

   Several subrutines were adapted from the C code used in
   the R package stats for the function phyper

 $Id: hyper.c,v 1.16 2009/01/09 15:18:46 neuhausi Exp $
**********************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define EXTERN
#include "hyper.h"

void
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **list_file,
		char **output_file,
		char **matrix_file,
		int *cols,
		int *rows,
		int *n_objects,
		char **obj_name,
		int *hits,
		double *Q,
		enum TAIL *tail)
{
  
  int c;
  extern char *optarg;
  int errflg = 0;
  
  *data_file = NULL;
  
  *list_file = NULL;

  *output_file = NULL;

  *matrix_file = NULL;

  *obj_name = NULL;
  
  *cols = 0;
  
  *rows = 0;
  
  *n_objects = 0;

  *hits = 0;

  *Q = 0.05;

  *tail = AUTO;

  while ((c = getopt(argc, argv, "d:l:o:m:c:r:k:n:i:q:t:vD")) != EOF)
    switch (c) {
    case 'd':
      *data_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*data_file, optarg);
      break;
    case 'l':
      *list_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*list_file, optarg);
      break;
    case 'o':
      *output_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*output_file, optarg);
      break;
    case 'm':
      *matrix_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*matrix_file, optarg);
      break;
    case 'c':
      *cols = atoi(optarg);
      break;
    case 'r':
      *rows = atoi(optarg);
      break;
    case 'k':
      *n_objects = atoi(optarg);
      break;
    case 'n':
      *obj_name = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*obj_name, optarg);
      break;
    case 'i':
      *hits = atoi(optarg);
      break;
    case 'q':
      *Q = atof(optarg);
      break;
    case 't':
      if (strcmp(optarg, "auto") == 0) {
	*tail = AUTO;
      } else if (strcmp(optarg, "lower") == 0) {
	*tail = LOWER;
      } else if (strcmp(optarg, "upper") == 0) {
	*tail = UPPER;
      } else {
	errflg++;
      }
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

  if (*data_file == NULL || (list_file = NULL && obj_name == NULL)) { 
    errflg++;
  }

  if (errflg) {
    printf("\nhyper: Hyper geometric distribution\n");
    printf("\nUsage: hyper <options> (options in square brackets are optional)\n\n");
    printf(" -d  <string>    data file. A tab delimeted file with the following format.\n");
    printf("                 A header line with variable names. Lines with ranked values\n");
    printf("                 preceded by the object name at the begining of the line.\n");
    printf(" -l  [<string>]  list file. The file containing the names of the variables. There\n");
    printf("                 must be one variable per line.\n");
    printf(" -o  [<string>]  output file. If no file is specified it defaults to stdout.\n");
    printf(" -m  [<string>]  matrix file. A file with pairwise distances between the top hits or\n");
    printf("                 an input file with list identifiers for pairwise comparisons.\n");
    printf(" -c  [<integer>] number of cols. It should correspond to the number of columns\n");
    printf("                 in the data file.\n");
    printf(" -r  [<integer>] number of rows. It should correspond to the number of rows\n");
    printf("                 in the data file.\n");
    printf(" -k  [<string>]  number of objects in the list file.\n");
    printf(" -n  [<string>]  object name. Run hyper for the object name.\n");
    printf(" -i  [<integer>] the number of hits to include in the matrix file.\n");
    printf(" -q  [<float>]   q value to use in false discovery rate. (default is 0.05).\n");
    printf(" -t  [<string>]  The tail to use for the p-value calculation. The available values\n");
    printf("                 are auto, lower or upper. The default is auto\n\n");
    printf("Examples\n");
    printf("\thyper -d ./t/hyper.small.dat -l ./t/hyper.lis\n");
    printf("\thyper -d ./t/hyper.big.dat -n GeneList2\n");
    printf("\thyper -d ./t/hyper.small.dat -m ./t/hyper.matrix.txt\n");
    printf("\n");
    exit(0);
  }
   
}

RESULT *
set_result(XDATA *xdata)
{

  int i, len;

  RESULT *res;
  VECTOR *vr;

  res = TYPE_ALLOC(RESULT);
  res->n = xdata->n;
  res->a = 7;
  res->attributes = ARRAY_ALLOC(res->a, char *);
  res->variables = ARRAY_ALLOC(res->n, char *);
  res->vectors = ARRAY_ALLOC(res->n, VECTOR *);

  // Allocate space for the actual data results
  for (i=0;i<res->n;i++) {
    vr = res->vectors[i] = TYPE_ALLOC(VECTOR);
    vr->data = ARRAY_ALLOC(res->a, double);
    vr->n = res->a;
  }

  // Set the attribute names
  len = strlen("q");
  res->attributes[0] = ARRAY_ALLOC((len+1), char);  
  strncpy(res->attributes[0], "q", len);
  res->attributes[0][len] = '\0';

  len = strlen("m");
  res->attributes[1] = ARRAY_ALLOC((len+1), char);  
  strncpy(res->attributes[1], "m", len);
  res->attributes[1][len] = '\0';

  len = strlen("n");
  res->attributes[2] = ARRAY_ALLOC((len+1), char);  
  strncpy(res->attributes[2], "n", len);
  res->attributes[2][len] = '\0';

  len = strlen("k");
  res->attributes[3] = ARRAY_ALLOC((len+1), char);  
  strncpy(res->attributes[3], "k", len);
  res->attributes[3][len] = '\0';

  len = strlen("odds ratio");
  res->attributes[4] = ARRAY_ALLOC((len+1), char);  
  strncpy(res->attributes[4], "odds ratio", len);
  res->attributes[4][len] = '\0';

  len = strlen("p-value");
  res->attributes[5] = ARRAY_ALLOC((len+1), char);  
  strncpy(res->attributes[5], "p-value", len);
  res->attributes[5][len] = '\0';

  len = strlen("fdr");
  res->attributes[6] = ARRAY_ALLOC((len+1), char);  
  strncpy(res->attributes[6], "fdr", len);
  res->attributes[6][len] = '\0';

  return res;

}

int
search_object(XDATA *xdata,
	      char *obj_name)
{

  int i, len;

  char c_obj_name[1000];
  char var_name[1000];

  // We compare everything in upper case

  len = strlen(obj_name);
  strncpy(c_obj_name, obj_name, len);
  c_obj_name[len] = '\0';
  uppercase(c_obj_name, len);
  
  for (i = 0; i < xdata->n; i++) {
    len = strlen(xdata->variables[i]);
    strncpy(var_name, xdata->variables[i], len);
    var_name[len] = '\0';
    uppercase(var_name, len);
    if (strcmp(c_obj_name, var_name) == 0) {
      return i;
    }
  }
  
  return (int) -1;

}

void
uppercase(char *str,
	  int len)
{

  int i;

  for (i = 0; i < len; i++) {
    str[i] = toupper((unsigned char)str[i]);
  }

}

void
set_vector(XDATA *xdata,
	   char **objects,
	   int n_objects,
	   VECTOR *v)
{ 

  int i, j, len;
  int idxs[n_objects];

  char sample_name[1000];
  char *uobjects[n_objects];

  for (i = 0; i < n_objects; i++) {
    idxs[i] = -1;
  }

  // We compare everything in upper case
  
  for (i = 0; i < n_objects; i++) {
    len = strlen(objects[i]);
    uobjects[i] = ARRAY_ALLOC((len+1), char);
    strncpy(uobjects[i], objects[i], len);
    uobjects[i][len] = '\0';
    uppercase(uobjects[i], len);
  }

  for (i = 0; i < xdata->s; i++) {
    len = strlen(xdata->samples[i]);
    strncpy(sample_name, xdata->samples[i], len);
    sample_name[len] = '\0';
    uppercase(sample_name, len);
    for (j = 0; j < n_objects; j++) {
      if (strcmp(uobjects[j], sample_name) == 0) {
	idxs[j] = i;
      }
    }
    if (idxs[j] < 0 && objects[j]) {
      fprintf(stderr, "Unable to find %s in data\n", objects[j]);
      exit(1);
    }
  }

  for (i = 0; i < xdata->s; i++) {
    v->data[i] = 0;
  }
  
  for (i = 0; i < n_objects; i++) {
    v->data[idxs[i]] = 1;  
  }

}

RESULT *
calculate_hyper (XDATA *xdata,
		 char *obj_name,
		 char **objects,
		 int n_objects,
		 int *hits,
		 double Q,
		 enum TAIL tail)
{

  int i, j, l, len, nsig;
  double q, m, n, k, e;
  double *data, *qdata;  
  double cutoff;
  RESULT *res;
  VECTOR *v, *v1, *vr;

  v = TYPE_ALLOC(VECTOR);
  v->n = 1;
  v->data = ARRAY_ALLOC(xdata->s, double);

  res = set_result(xdata);

  data = ARRAY_ALLOC(res->n, double);
  qdata = ARRAY_ALLOC(res->n, double);

  // We do some validation and load the
  // data to use as query into a new vector
  if (obj_name != NULL) {
    l = search_object(xdata, obj_name);
    if (l < 0) {
      fprintf(stderr, "Unable to find %s in data\n", obj_name);
      exit(1);
    }
    v1 = xdata->vectors[l];
    for (i = 0; i < xdata->s; i++) {
      v->data[i] = v1->data[i];  
    }
  } else {
    set_vector(xdata, objects, n_objects, v);
  }

  k = 0.0;
  for (i = 0; i < xdata->s; i++) {
    if (v->data[i] > 0) {
      k++;
    }
  }

  // We compute the hyper-geometric dist
  for (i = 0; i < xdata->n; i++) {
    //Get the variable name
    len = strlen(xdata->variables[i]);
    res->variables[i] = ARRAY_ALLOC((len+1), char);
    strncpy(res->variables[i], xdata->variables[i], len);
    res->variables[i][len] = '\0';
    // Get the data in each vector
    v1 = xdata->vectors[i];
    q = 0.0;
    m = 0.0;
    e = 0.0;
    for (j = 0; j < xdata->s; j++) {
      if (v->data[j] > 0 && v1->data[j] > 0) {
	q++;
      }
      if (v->data[j] < 0 && v1->data[j] != 1) {
	e++;
      }
      if (v1->data[j] > 0) {
	m++;
      }
    }
    n = xdata->s - (m + e); 
    hyper(q, m, n, k, res, i, tail);
    vr = res->vectors[i]; 
    data[i] = vr->data[5];    
  }

  // Calculate the fdr
  fdr(res->n, data, qdata, Q, 1, &cutoff, &nsig);

  // Load the data into the result vector
  for (i = 0; i < res->n; i++) {
    vr = res->vectors[i]; 
    vr->data[6] = qdata[i];
  }

  // Modifiy the hits
  if (*hits == 0) {
    *hits = nsig;
  }

  free(data);
  free(qdata);
  free(v);

  return res;

}

void
hyper (double q,
       double m,
       double n,
       double k,
       RESULT *res,
       int i,
       enum TAIL tail)
{

  double oldn; 
  double d;
  double pd;
  double pval;
  VECTOR *vr;

  vr = res->vectors[i];

  vr->data[0] = q;
  vr->data[1] = m;
  vr->data[2] = n;
  vr->data[3] = k;

  if ((m-q) == 0 || (k-q) == 0) {
    vr->data[4] = 999999; 
  } else {
    vr->data[4] = (q*(n-k+q))/((m-q)*(k-q));
  }

  if ((q*(m+n)) > (k*m)) {
    oldn = n;
    n = m;
    m = oldn;
    q = k-q; 
    if (tail != LOWER) {
      tail = LOWER;
    } else if (tail == LOWER) {
      tail = UPPER;
    }
  } else {
    if (tail == AUTO) {
      tail = LOWER;
    };
  }

  if (tail == LOWER) {
    d  = dhyper(q, m, n, k);
    pd = pdhyper(q, m, n, k);
    pval = d * pd;
  } else if (tail == UPPER) {
    d  = dhyper(q-1, m, n, k);
    pd = pdhyper(q-1, m, n, k);
    pval = 0.5 - d * pd + 0.5;
  }

  vr->data[5] = pval;
  
}

double
dhyper (double q,
	double m,
	double n,
	double k) {

  double p, pq, p1, p2, p3;

  p = ((double)k)/((double)(m+n));
  pq = ((double)(m+n-k))/((double)(m+n));

  p1 = dbinom_raw(q,  m, p,pq);
  p2 = dbinom_raw(k-q,n, p,pq);
  p3 = dbinom_raw(k,m+n, p,pq);

  return(p1*p2/p3);

}


double
dbinom_raw (double q, 
            double k, 
            double p, 
            double pq) {

  double lf, lc;
 
  if (p == 0) return ((q == 0) ? 1 : 0);
  if (pq == 0) return ((q == k) ? 1 : 0);
  
  if (q == 0) {
    if (k == 0) return 1;
    lc = (p < 0.1) ? -bd0(k,k*pq) - k*p : k*log(pq);
    return(exp(lc));
  }
  if (q == k) {
    lc = (pq < 0.1) ? -bd0(k,k*p) - k*pq : k*log(p);
    return(exp(lc));
  }
  if (q < 0 || q > k) return(0);
                                                                                                                                            
  lc = stirlerr(k) - stirlerr(q) - stirlerr(k-q) - bd0(q,k*p) - bd0(k-q,k*pq);
  
  lf = log(PIX2) + log(q) + log(k-q) - log(k);
  
  return exp(lc - 0.5*lf);
}
 

double 
stirlerr (double n) {
  
  double nn;

  #define S0 0.083333333333333333333       /* 1/12 */
  #define S1 0.00277777777777777777778     /* 1/360 */
  #define S2 0.00079365079365079365079365  /* 1/1260 */
  #define S3 0.000595238095238095238095238 /* 1/1680 */
  #define S4 0.0008417508417508417508417508/* 1/1188 */
 
  const double sferr_halves[31] = {
    0.0, /* n=0 - wrong, place holder only */
    0.1534264097200273452913848,  /* 0.5 */
    0.0810614667953272582196702,  /* 1.0 */
    0.0548141210519176538961390,  /* 1.5 */
    0.0413406959554092940938221,  /* 2.0 */
    0.03316287351993628748511048, /* 2.5 */
    0.02767792568499833914878929, /* 3.0 */
    0.02374616365629749597132920, /* 3.5 */
    0.02079067210376509311152277, /* 4.0 */
    0.01848845053267318523077934, /* 4.5 */
    0.01664469118982119216319487, /* 5.0 */
    0.01513497322191737887351255, /* 5.5 */
    0.01387612882307074799874573, /* 6.0 */
    0.01281046524292022692424986, /* 6.5 */
    0.01189670994589177009505572, /* 7.0 */
    0.01110455975820691732662991, /* 7.5 */
    0.010411265261972096497478567, /* 8.0 */
    0.009799416126158803298389475, /* 8.5 */
    0.009255462182712732917728637, /* 9.0 */
    0.008768700134139385462952823, /* 9.5 */
    0.008330563433362871256469318, /* 10.0 */
    0.007934114564314020547248100, /* 10.5 */
    0.007573675487951840794972024, /* 11.0 */
    0.007244554301320383179543912, /* 11.5 */
    0.006942840107209529865664152, /* 12.0 */
    0.006665247032707682442354394, /* 12.5 */
    0.006408994188004207068439631, /* 13.0 */
    0.006171712263039457647532867, /* 13.5 */
    0.005951370112758847735624416, /* 14.0 */
    0.005746216513010115682023589, /* 14.5 */
    0.005554733551962801371038690  /* 15.0 */
  };

 
  if (n <= 15.0) {
    nn = n + n;
    if (nn == (int)nn) {
      return (sferr_halves[(int)nn]);
    } else {
      fprintf(stderr, "%f not integer\n", nn);
      exit(1);
    }
  } else {
    nn = n*n;
    if (n>500) return ((S0-S1/nn)/n);
    if (n> 80) return ((S0-(S1-S2/nn)/nn)/n);
    if (n> 35) return ((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
    /* 15 < n <= 35 : */
    return ((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);
  }

}

double 
bd0 (double x, 
     double np) {

  double ej, s, s1, v;
  int j;
  
  if ((abs(x-np)) < (0.1*(x+np))) {
    v = (x-np)/(x+np);
    s = (x-np)*v; /* s using v -- change by MM */
    ej = 2*x*v;
    v = v*v;
    for (j=1; ; j++) { /* Taylor series */
      ej *= v;
      s1 = s+ej/((j<<1)+1);
      if (s1==s) { /* last term was effectively 0 */
	return (s1);
      }
      s = s1;
    }
  } else {
    // | x - np |  is not too small */
    return (x*log(x/np)+np-x);
  }

}                                                                                                                                                                                                     

double
pdhyper (double q,
	 double m,
	 double n,
	 double k) {
    
  double sum = 0;
  double term = 1;
 
  while (q > 0 && term >= DBL_EPSILON * sum) {
    term *= q * (n - k + q) / (k + 1 - q) / (m + 1 - q);
    sum += term;
    q--;
  }
 
  return 1 + sum;

}

RESULT *
set_result_matrix(XDATA *xdata,
		  int *idx,
		  int hits)
{

  int i, len;
  RESULT *resmat;
  VECTOR *vr;

  resmat = TYPE_ALLOC(RESULT);
  resmat->n = hits;
  resmat->a = hits;
  resmat->attributes = ARRAY_ALLOC(resmat->a, char *);
  resmat->variables = ARRAY_ALLOC(resmat->n, char *);
  resmat->vectors = ARRAY_ALLOC(resmat->n, VECTOR *);

  // Allocate space for the actual data results
  for (i = 0; i < resmat->n; i++) {
    vr = resmat->vectors[i] = TYPE_ALLOC(VECTOR);
    vr->data = ARRAY_ALLOC(resmat->a, double);
    vr->n = resmat->a;
  }

  // Set up the row and column labels
  for (i = 0; i < hits; i++) {
    len = strlen(xdata->variables[idx[i]]);
    resmat->attributes[i] = ARRAY_ALLOC((len+1), char);  
    resmat->variables[i] = ARRAY_ALLOC((len+1), char);  
    strncpy(resmat->attributes[i],  xdata->variables[idx[i]], len);
    strncpy(resmat->variables[i],  xdata->variables[idx[i]], len);
    resmat->attributes[i][len] = '\0';
    resmat->variables[i][len] = '\0';
  }

  return resmat;

}

void
calculate_pw_hyper(XDATA *xdata,
		   RESULT *res,
		   int hits,
		   char **ids,
		   double Q,
		   enum TAIL tail,
		   RESULT **resmatp,
		   RESULT **resmato,
		   int exist)
{

  int i, j, ii, jj, s, perm[res->n];
  double q, m, n, k, e;
  double tmp[res->n];
  VECTOR *v, *tmpv, *v1, *v2, *vr, *vmatp, *vmato, *vmatp1, *vmato1;
  RESULT *tmpres;

  // Temporary reult object
  tmpres = TYPE_ALLOC(RESULT);
  tmpres->n = 1;
  tmpres->a = 7;
  tmpres->vectors = ARRAY_ALLOC(res->n, VECTOR *);
  tmpv = tmpres->vectors[0] = TYPE_ALLOC(VECTOR);
  tmpv->data = ARRAY_ALLOC(tmpres->a, double);
  tmpv->n = tmpres->a;

  if (exist > 0) {
    find_idx_for_ids(xdata, ids, hits, perm);
  } else {
  // Sort the results by p-value
    for (i = 0; i < res->n; i++) {
      v = res->vectors[i];
      tmp[i] = v->data[5];		
    }
    sort_by_index(res->n, tmp, perm);
  }

  *resmatp = set_result_matrix(xdata, perm, hits);
  *resmato = set_result_matrix(xdata, perm, hits);

  // Calculate the pairwise hyper for the top triangle
  for (i = 0; i < hits; i++) {
    // printf("Current iteration: %d \n",i); 
    // vmat is where the results will be stored
    vmatp = (*resmatp)->vectors[i];
    vmato = (*resmato)->vectors[i];
    ii = perm[i];
    v1 = xdata->vectors[ii];
    k = 0.0;
    for (s = 0; s < xdata->s; s++) {
      if (v1->data[s] > 0) {
	k++;
      }
    }
    for (j = 0; j < hits; j++) {
      if (j >= i) {
	jj = perm[j];
	v2 = xdata->vectors[jj];
	q = 0.0;
	m = 0.0;
	e = 0.0;
	for (s = 0; s < xdata->s; s++) {
	  if (v1->data[s] > 0 && v2->data[s] > 0) {
	    q++;
	  }
	  if (v1->data[s] < 0 && v2->data[s] != 1) {
	    e++;
	  }
	  if (v2->data[s] > 0) {
	    m++;
	  }
	}
	n = xdata->s - (m + e); 
	hyper(q, m, n, k, tmpres, 0, tail);
	vr = tmpres->vectors[0];
	vmatp->data[j] = vr->data[5];
	vmato->data[j] = vr->data[4];
      }
    }
  }

  // Copy the reciprocal values to the bottom triangle
  for (i = 0; i < hits; i++) {
    vmatp = (*resmatp)->vectors[i];
    vmato = (*resmato)->vectors[i];
    for (j = hits - 1; j >= 0; j--) {
      vmatp1 = (*resmatp)->vectors[j];
      vmato1 = (*resmato)->vectors[j];
      if (j >= i) {
	vmatp1->data[i] = vmatp->data[j]; 
	vmato1->data[i] = vmato->data[j]; 
      }
    }
  }

}

void
find_idx_for_ids(XDATA *xdata,
		 char **ids,
		 int hits, 
		 int *perm)
{

  int i, j, len;

  char var_name[1000];
  char *uids[hits];

  // We compare everything in upper case

  for (i = 0; i < hits; i++) {
    len = strlen(ids[i]);
    uids[i] = ARRAY_ALLOC((len+1), char);
    strncpy(uids[i], ids[i], len);
    uids[i][len] = '\0';
    uppercase(uids[i], len);
  }

  for (i = 0; i < xdata->n; i++) {
    len = strlen(xdata->variables[i]);
    strncpy(var_name, xdata->variables[i], len);
    var_name[len] = '\0';
    uppercase(var_name, len);
    for (j = 0; j < hits; j++) {
      if (strcmp(uids[j], var_name) == 0) {
	perm[j] = i;
      }
    }
  }

  free(uids);

}

void
print_matrix_results(RESULT *resp,
		     RESULT *reso,
		     char *out_file)
{

  char *name;

  name = ARRAY_ALLOC((strlen(out_file) + strlen(".xxx") + 1), char);
  strcpy(name, out_file);
  strcat(name, ".pvl");
  print_results(name, resp, 0, 0);

  name = ARRAY_ALLOC((strlen(out_file) + strlen(".xxx") + 1), char);
  strcpy(name, out_file);
  strcat(name, ".odd");
  print_results(name, reso, 0, 0);

  free(name);

}

int 
main(int argc,
     char **argv)
{

  char *data_file, *list_file, *matrix_file, *output_file;
  int cols, rows, n_objects, hits, exist;
  enum TAIL tail;
  double Q;
  char *obj_name;
  char **objects;
  char **ids;
  XDATA *xdata;
  RESULT *res, *resmatp, *resmato;
  
  verbose = 0;

  debug = 0;

  parse_arguments(argc, argv,
		  &data_file,
		  &list_file,
		  &output_file,
		  &matrix_file,		  
		  &cols, &rows, &n_objects,
		  &obj_name, &hits, &Q, &tail);
    
  if (cols == 0) {
    cols = count_fields(data_file);
  }

  if (rows == 0) {
    rows = count_lines(data_file);
    rows--;
  }

  xdata = read_xdata(data_file, cols, rows);

  if (list_file != NULL) {
    if (n_objects == 0) {
      n_objects = count_lines(list_file);
    }
    objects = read_file_column(list_file, n_objects, 0);
  } else {
    n_objects = 0;
  }

  res = calculate_hyper(xdata, obj_name, objects, n_objects, &hits, Q, tail);

  if (hits > res->n) {
    hits = res->n;
  }

  print_results(output_file, res, 1, 5);
 
  if (matrix_file != NULL) {
    exist = file_exists(matrix_file);
    if (exist > 0) {
      hits = count_lines(matrix_file);
      ids = read_file_column(matrix_file, hits, 0);
    }
    calculate_pw_hyper(xdata, res,  hits, ids, Q, tail, &resmatp, &resmato, exist);    
    print_matrix_results(resmatp, resmato, matrix_file);
    free(resmatp);
    free(resmato);
    free(matrix_file);
  }

  if (n_objects > 0) {
    free(objects);
  }

  free(data_file);
  free(list_file);
  if (output_file != NULL) {
    free(output_file);
  }
  free(obj_name);
  free(xdata);
  free(res);

  exit(0);

}
