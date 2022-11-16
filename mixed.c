/*********************************************************************
 Copyright 2004 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

 $Id: mixed.c,v 1.44 2010/04/15 17:32:21 neuhausi Exp $
**********************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define EXTERN
#include "mixed.h"

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
		char ***fixed_effs,
		int *n_fixed_effs,
		char ***interacts,
		int *n_interacts,
		char ***random_facs,
		int *n_random_facs,
		enum METHOD *method,
		int *ssq,
		double *Q,
		char **reference,
		int *lsmo,
		int *islog,
		int *permuts)
{
  
  int c, arg_len, max_objs, max_vars;
  extern char *optarg;
  int errflg = 0;
  
  *data_file = NULL;
  
  *factor_file = NULL;
  
  *covariate_file = NULL;

  *output_file = NULL;

  *reference = NULL;

  *samples = 0;
  
  *variables = 0;
  
  /* Maximum of 100 objects */
  max_vars = 100;

  *variable_names = ARRAY_ALLOC(max_vars, char *);
  
  *variable_number = 0;

  *n_fixed_effs = 0;

  *n_interacts = 0;

  *n_random_facs = 0;

  *method = REML;

  *ssq = 3;

  max_objs = 10;

  *fixed_effs = ARRAY_ALLOC(max_objs, char *);
  
  *interacts = ARRAY_ALLOC(max_objs, char *);

  *random_facs = ARRAY_ALLOC(max_objs, char *);
  
  *Q = 0.05;

  *lsmo = 0;

  *islog = 0;

  *permuts = 10;
  
  while ((c = getopt(argc, argv, "d:f:c:o:s:g:n:m:i:r:e:t:p:q:x:lzvD")) != EOF)
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
      if (*n_fixed_effs <= max_objs) {
	arg_len = strlen(optarg);
	(*fixed_effs)[*n_fixed_effs] = ARRAY_ALLOC((arg_len+1), char);
	strcpy((*fixed_effs)[*n_fixed_effs], optarg);
	*n_fixed_effs += 1;
      }
      break;
    case 'i':
      if (*n_interacts <= max_objs) {
	arg_len = strlen(optarg);
	(*interacts)[*n_interacts] = ARRAY_ALLOC((arg_len+1), char);
	strcpy((*interacts)[*n_interacts], optarg);
	*n_interacts += 1;
      }
      break;
    case 'r':
      if (*n_random_facs <= max_objs) {
	arg_len = strlen(optarg);
	(*random_facs)[*n_random_facs] = ARRAY_ALLOC((arg_len+1), char);
	strcpy((*random_facs)[*n_random_facs], optarg);
	*n_random_facs += 1;
      }
      break;
    case 'e':
      if (strcmp(optarg, "reml") == 0) {
	*method = REML;
      } else if (strcmp(optarg, "ml") == 0) {
	*method = ML;
      } else if (strcmp(optarg, "mivque") == 0) {
	*method = MIVQUE;
      }
      break;
    case 't':
      *ssq = atoi(optarg);
      if (*ssq < 2 || *ssq > 3) {
	errflg++;
      }
      break;
    case 'q':
      *Q = atof(optarg);
      break;
    case 'x':
      *reference = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*reference, optarg);
      break;
    case 'l':
      *lsmo = 1;
      break;
    case 'z':
      *islog = 1;
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

  if (*data_file == NULL || (*factor_file == NULL && *covariate_file == NULL)) { 
    errflg++;
  }
  if (*n_fixed_effs+*n_interacts==0) {
    errflg++;
  }

  if (errflg) {
    printf("\nGeneral Linear Model / Mixed Linear Model: A simple program to run AN[C]OVAs.");
    printf("Supports all sort of experimental designs (at least those that I can imagine, balanced,\n");
    printf("unbalanced, nested, empty cells, missing treatment combinations, etc). The program\n");
    printf("automatically detects nested factors so there is no need to explicitly specify them in\n");
    printf("the model. It can handle up to nine factors and nine covariates.\n");
    printf("\nUsage: glm <options> (options in square brackets are optional)\n\n");
    printf(" -d  <string>    data file. A tab delimeted file with the following format:\n");
    printf("                 A header line with the sample names followed by one or more lines\n");
    printf("                 with the data. The first column of these lines must start with the\n");
    printf("                 descriptor name. Missing data must be entered as NA.\n");
    printf(" -f  <string>    factor file. A tab delimited file containing the classes for\n");
    printf("                 the samples.\n");
    printf(" -c  [<string>]  covariate file. A tab delimited file containing the covariates for\n");
    printf("                 the samples.\n");
    printf(" -o  [<string>]  output file. If no file is specified it defaults to stdout.\n");
    printf(" -s  [<integer>] number of samples.\n");
    printf(" -g  [<integer>] number of variables.\n");
    printf(" -n  <string>    variable name. Run glm only for the specified variable name\n");
    printf(" -m  <string>    fixed effect(s) name (must match one name in the factor file).\n");
    printf("                 Multiple fixed effects are entered with multiple -m in the command\n");
    printf("                 line. Do NOT enter interactions here.\n");
    printf(" -i  [<string>]  interaction(s) between factors (must match the names in the factor\n");
    printf("                 file). Multiple interactions are entered with multiple -i in the command\n");
    printf("                 line. Each interaction must be joined by an asterix '*', for example,\n");
    printf("                 subject*treatment.\n");
    printf(" -e  [<string>]  [reml|ml|mivque] method to calculate mixed models (defaults to reml).\n");
    printf(" -r  [<string>]  random factor name (must match one of the names in the factor file).\n");
    printf("                 Multiple random factors are entered with multiple -r in the command\n");
    printf("                 line.\n");
    printf(" -l  <switch>    output least square means\n");
    printf(" -z  <switch>    data is log transformed (only used when calculating ratios/log ratios\n");
    printf(" -x  <string>    contrast to calculate. It needs to identify a factor and a level\n");
    printf("                 separated by a colon ':'. An example would be 'Treatment:Control',\n");
    printf("                 where Treatment is a factor and Control is a level of Treatment.\n");
    printf(" -t  [<integer>] type of sum of squares [2|3]. Type III is implemented only\n");
    printf(" -q  [<float>]   q value to use in false discovery rate. (default is 0.05).\n\n");
    printf(" -p  [<integer>] number of random permutations (not implemented yet).\n\n");
    printf("Examples\n");
    printf("\tmixed -d ./t/mixed.dat -f ./t/mixed.fac -m A -m B -r B -v -n var\n");
    printf("\n");
    exit(0);
  }
   
}

RESULT *
set_result(int ndata,
	   FACTOR *factor,
	   COVARIATE *covariate,
	   MODEL *model,
	   MATRIX *d,
	   double Q,
	   char *reference,
	   int lsmo)
{

  int i, j, f, t, r, n, idx, len;
  char str[1000], stra[1000], strn[100];
  RESULT *res;
  VECTOR *vr;

  res = TYPE_ALLOC(RESULT);
  res->n = ndata;
  if (model->r > 0) {
      r = 0;
      for (i=0;i<model->i;i++) {
	  if (model->randinter[i] > 0) {
	      r++;
	  }
      }
      res->a = (2 * model->f) + (2 * model->i) - (model->r + r) + 3;
  } else {
      res->a = model->f + model->i + 5;
  }
  if (lsmo > 0) {
    for (i=0;i<model->t;i++) {
      for (j=0;j<model->termlen[i];j++) {
        if (model->zeros[i][j] > -1) {
	  res->a++;
	}
      }
    }
  }
  if (reference != NULL) {
    for (i=0;i<model->t;i++) {
      for (j=0;j<model->termlen[i];j++) {
        if (model->zeros[i][j] > -1) {
	  res->a++;
        }
      }
    }
    res->a--;
  }
  res->attributes = ARRAY_ALLOC(res->a, char *);
  res->variables = ARRAY_ALLOC(res->n, char *);
  res->vectors = ARRAY_ALLOC(res->n, VECTOR *);
  for (i=0;i<ndata;i++) {
    vr = res->vectors[i] = TYPE_ALLOC(VECTOR);
    vr->data = ARRAY_ALLOC(res->a, double);
    vr->n = res->a;
  }

  idx=0;
  // Fixed Effects
  for (i=0;i<model->f;i++) {
      f = model->fixed[i];
      t = model->tyfixed[i];
      if (t==0) {
	  // factor
	  strcpy(str, factor->fnames[f]);
      } else {
	  // covariate
	  strcpy(str, covariate->cnames[f]);
      }
      if (model->r > 0) {
	  if (model->infactorrand[f] == 0) {
	      set_attribute_result(idx,str,res);
	      idx++;
	  }
      } else {
	  set_attribute_result(idx,str,res);
	  idx++;
      }
  }
  // Interactions
  for (i=0;i<model->i;i++) {
      n=0;
      strcpy(str, "");
      for (j=0;j<model->numinter[i];j++) {
	  f = model->inter[i][j];
	  t = model->tyinter[i][j];
	  if (t==0 && model->infactornest[f]>0) {
	      strcpy(strn, factor->fnames[f]);
	      n++;
	  } else {
	      if (t==1) {
		  // factor
		  strcat(str, covariate->cnames[f]); 
	      } else {
		  // covariate
		  strcat(str, factor->fnames[f]);
	      }
	      strcat(str, "*");
	  }
      }
      len = strlen(str)-1;
      strncpy(stra, str, len);
      stra[len]='\0';
      if (n>0) {
	  strcat(strn, "(");
	  strcat(strn, stra);
	  strcat(strn, ")");
	  if (model->randinter[i] == 0) {
	      set_attribute_result(idx,strn,res);
	      idx++;
	  }
      } else {
	  if (model->randinter[i] == 0) {
	      set_attribute_result(idx,stra,res);
	      idx++;
	  }
      }
  }
  // Variance Components from Random factors
  if (model->r > 0) {
      // Random effects
      for (i=0;i<model->f;i++) {
	  f = model->fixed[i];
	  t = model->tyfixed[i];
	  if (model->infactorrand[f]>0) {
	      strcpy(str, "Var(");
	      if (t==0) {
		  // factor
		  strcat(str, factor->fnames[f]);
	      } else {
		  // covariate
		  strcat(str, covariate->cnames[f]);
	      }
	      strcat(str, ")");
	      set_attribute_result(idx,str,res);
	      idx++;
	  }
      }
      // Random Interactions
      for (i=0;i<model->i;i++) {
	  if (model->randinter[i] > 0) {
	      n=0;
	      strcpy(str, "");
	      for (j=0;j<model->numinter[i];j++) {
		  f = model->inter[i][j];
		  t = model->tyinter[i][j];
		  if (t==0 && model->infactornest[f]>0) {
		      strcat(strn, factor->fnames[f]);
		      n++;
		  } else {
		      if (t==1) {
			  // factor
			  strcat(str, covariate->cnames[f]); 
		      } else {
			  // covariate
			  strcat(str, factor->fnames[f]);
		      }
		      strcat(str, "*");
		  }
	      }
	      len = strlen(str)-1;
	      strncpy(stra, str, len);
	      stra[len]='\0';
	      strcpy(str, "Var(");
	      if (n>0) {
		  strcat(strn, "(");
		  strcat(strn, stra);
		  strcat(strn, ")");
		  strcat(str, strn);
	      } else {
		  strcat(str, stra);
	      }
	      strcat(str, ")");
	      set_attribute_result(idx,str,res);
	      idx++;
	  }
      }
  }

  if (model->r > 0) {
      strcpy(str, "Var(Error)");
      set_attribute_result(idx,str,res);
      idx++;
      strcpy(str, "AIC");
      set_attribute_result(idx,str,res);
      idx++;
      strcpy(str, "-2LogLik");
      set_attribute_result(idx,str,res);
      idx++;
      for (i=0;i<model->f;i++) {
	  if (model->infactorrand[i] == 0) {
	      f = model->fixed[i];
	      t = model->tyfixed[i];
	      strcpy(str, "FDR (");
	      sprintf(strn, "%5.3f", Q);
	      strcat(str, strn);
	      strcat(str, ") ");
	      if (t==0) {
		  // factor
		  strcat(str, factor->fnames[f]);
	      } else {
		  // covariate
		  strcat(str, covariate->cnames[f]);
	      }
	      set_attribute_result(idx,str,res);
	      idx++;
	  }
      }
      for (i=0;i<model->i;i++) {
	  if (model->randinter[i] == 0) {
	      n=0;
	      strcpy(str, "FDR (");
	      sprintf(strn, "%5.3f", Q);
	      strcat(str, strn);
	      strcat(str, ") ");
	      for (j=0;j<model->numinter[i];j++) {
		  f = model->inter[i][j];
		  t = model->tyinter[i][j];
		  if (t==0 && model->infactornest[f]>0) {
		      strcpy(strn, factor->fnames[f]);
		      n++;
		  } else {
		      if (t==1) {
			  // factor
			  strcat(str, covariate->cnames[f]); 
		      } else {
			  // covariate
			  strcat(str, factor->fnames[f]);
		      }
		      strcat(str, "*");
		  }
	      }
	      len = strlen(str)-1;
	      strncpy(stra, str, len);
	      stra[len]='\0';
	      if (n>0) {
		  strcat(strn, "(");
		  strcat(strn, stra);
		  strcat(strn, ")");
		  set_attribute_result(idx,strn,res);
		  idx++;
	      } else {
		  set_attribute_result(idx,stra,res);
		  idx++;
	      }
	  }
      }
  } else {
      strcpy(str, "Model");
      set_attribute_result(idx,str,res);
      idx++;
      strcpy(str, "R-square");
      set_attribute_result(idx,str,res);
      idx++;
      strcpy(str, "Adj R-sq");
      set_attribute_result(idx,str,res);
      idx++;
      strcpy(str, "AIC");
      set_attribute_result(idx,str,res);
      idx++;
      strcpy(str, "FDR (");
      sprintf(strn, "%5.3f", Q);
      strcat(str, strn);
      strcat(str, ")");
      set_attribute_result(idx,str,res);
      idx++;
  }

  if (lsmo > 0) {
    for (i=0;i<model->t;i++) {
      for (j=0;j<model->termlen[i];j++) {
        if (model->zeros[i][j] > -1) {
	  strcpy(str, "LSM ");
	  strcat(str, model->factors[i]);
	  strcat(str, ":");
	  strcat(str, model->levels[i][j]);
	  len = strlen(str);
	  str[len]='\0';
	  set_attribute_result(idx,str,res);
	  idx++;
	}
      }
    }
  }

  if (reference != NULL) {
    for (i=0;i<model->t;i++) {
      for (j=0;j<model->termlen[i];j++) {
        if (model->zeros[i][j] > -1) {
	  strcpy(str, model->factors[i]);
	  strcat(str, ":");
	  strcat(str, model->levels[i][j]);
	  len = strlen(str);
	  str[len]='\0';
	  if (strcmp(str, reference)) {
	    strcpy(stra, str);
	    strcat(stra, " / ");
	    strcat(stra, reference);
	    set_attribute_result(idx,stra,res);
	    idx++;
	  }
	}
      }
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

void
validate_reference(MODEL *model,
		   char *reference)
{

  int i, j, len;
  char str[1000];

  for (i=0;i<model->t;i++) {
    for (j=0;j<model->termlen[i];j++) {
      strcpy(str, model->factors[i]);
      strcat(str, ":");
      strcat(str, model->levels[i][j]);
      len = strlen(str);
      str[len]='\0';
      if (!strcmp(str, reference)) {
	return;
      }
    }
  }

  fprintf(stderr, "Contrast %s not found in model\n", reference);
  exit(1);

} 

void
mult_test_correction(RESULT *res,
		     MODEL *model,
		     double Q)
{

    int i, j, s, e, nsig;
    VECTOR *v;
    double data[res->n], qdata[res->n];  
    double cutoff;

  if (model->r > 0) {
      s = model->t+3;
      e = s + model->f - model->r;
      for (i=0;i<model->i;i++) {
	  if (model->randinter[i] == 0) {
	      e++;
	  }
      }
      for (i=s;i<e;i++) {
	  for (j=0;j<res->n;j++) {
	      v = res->vectors[j];
	      data[j] = v->data[i-s];

	  }
	  fdr(res->n,data,qdata,Q,1,&cutoff,&nsig);
	  for (j=0;j<res->n;j++) {
	      v = res->vectors[j];
	      v->data[i] = qdata[j];
	  }
      }
  } else {
      for (i=0;i<res->n;i++) {
	  v = res->vectors[i];
	  data[i] = v->data[model->t];
      }
      fdr(res->n,data,qdata,Q,1,&cutoff,&nsig);
      for (i=0;i<res->n;i++) {
	  v = res->vectors[i];
	  v->data[model->f + model->i + 4] = qdata[i];
      }
  }

}

MATRIX *
set_G_derivatives(MATRIX *G,
		  int **indeces,
		  int index)
{

  int i;
  MATRIX *deriv;

  deriv = fill(G->n, G->m, 0.0);
  for (i = indeces[index][0]; i < indeces[index][1]; i++) {
    deriv->data[i][i] = 1;
  }

  return deriv;

}

void
update_G_matrix(MATRIX **G,
	     MATRIX *parameters,
		int **indeces)
{

  // Populate the G matrix awith the estimated
  // parameters from REML and ML or MIVQUE

  int i, j;

  for (i = 0; i < parameters->n - 1; i++) {
    for (j = indeces[i][0]; j < indeces[i][1]; j++) {
      (*G)->data[j][j] = parameters->data[i][0];
    }
  }

  //(*G) = multiply_number(K,parameters->data[0][0],0);
  
}

void
update_R_matrix(MATRIX **R,
		MATRIX *parameters)
{

  int i;

  // Populate the R matrix awith the estimated
  // parameters from REML and ML or MIVQUE

  for (i = 0; i < (*R)->n; i++) {
    (*R)->data[i][i] = parameters->data[parameters->n - 1][0];
  }

}

RESULT *
calculate_model(XDATA *xdata,
		FACTOR *factor,
		COVARIATE *covariate,
		char **variable_names,
		int variable_number,
		MODEL *model,
		enum METHOD method,
		int ssq,
		double Q,
		char *reference,
		int lsmo,
		int islog,
		int permuts)
{

  int i, ii, j, k, l, len, n, g, s, na, nrand;
  int *variables, **indeces, **findeces, **rindeces, *rank, *rna, *ndf;
  RESULT *res;
  MATRIX *des, *desf, *desr, *gef, *X, *R, *G;
  MATRIX *desna, *desfna, *desrna, *gefna, *Xna,*Rna;
  MATRIX *hess, *grad;
  NMATRIX *G_deriv, *R_deriv;
  NMATRIX *R_derivna;
  NMATRIX *hyp, *con, *lsm;
  NMATRIX *hypna, *conna, *lsmna;
  VECTOR *v, *v1;

  // Allocate memory for an array with the indeces for the data
  // vectors. If there is one variable the lenght will be one otherwise
  // the lenght will be the size of the number of variables in the xdata
  // object.
  if (variable_number > 0) {
    n = variable_number;
  } else {
    n = xdata->n;   
  }
  variables = ARRAY_ALLOC(n, int);  
      
  // The number of samples. If there are factors we get the length
  // from the factor object otherwise we get it from the covariate
  // object
  s = factor->s > 0 ? factor->s : covariate->s;

  // An array to store the data with NA values so we can clean the
  // design matrix if necesary
  rna = ARRAY_ALLOC(s, int);

  // Generate some matrices used in all all calculations to speed up
  // the glm implementation (for data without NA values)
  
  // This is the design matrix. We will always create it even if there
  // are random factors since we will use it to calculate LSMEANS
  des = set_design_matrix(factor, covariate, model, &indeces);

  // These are the numerator degrees of freedom
  ndf = calculate_dof(des, model, indeces, ssq);

  // If there are random factors we need to set an additional series
  // of matrices
  if (model->r > 0) {
  
    // Number of random terms including random interactions
    nrand = 0;
    for (i=0;i<model->t;i++) {
       if (model->randdep[i] > 0) {
         nrand++;
       }   
    }
 
    // A matrix with only the fixed effects
    desf = set_fixed_matrix(factor, covariate, model, &findeces);

    // A matrix with only the random effects
    desr = set_random_matrix(factor, covariate, model, &rindeces);

    //Gradient Matrix
    grad = fill(nrand+1, 1, 0);

    // Hessian Matrix
    hess = fill(nrand+1, nrand+1, 0);
    // Set up the R matrix which contains the error
    R = identity(model->s - model->k);

    // Set the G matrix that contains the random factors
    G = identity(desr->m);
 
    // Set the G Derivative matrices
    G_deriv = TYPE_ALLOC(NMATRIX);
    G_deriv->n =  nrand;
    G_deriv->m = ARRAY_ALLOC(G_deriv->n, MATRIX *);
    for (i = 0; i < G_deriv->n; i++) {
      G_deriv->m[i] = set_G_derivatives(G,rindeces,i);
    }
      
    // Set the R Derivative matrices
    R_deriv = TYPE_ALLOC(NMATRIX);
    R_deriv->n =  1;
    R_deriv->m = ARRAY_ALLOC(R_deriv->n, MATRIX *);
    for (i = 0; i < 1; i++) {
      R_deriv->m[i] = identity(model->s - model->k);
    }
    // Set the Estimable function matrix
    gef = set_gef_matrix(desf, &rank);

    // Type of sum of squares
    if (ssq == 2) {
      hyp = set_ef3_matrix(model, gef, rank, findeces, 0);
    } else {
      hyp = set_ef3_matrix(model, gef, rank, findeces, 0);
    }

    // The Least Square Means Matrix
      
    lsm = set_lsmeans_matrix(model, des, indeces);

    // inv(X'X) * X'
    X = set_X_matrix(desf);
 
    // C * inv(X'X) * C'
    
    con = set_contrast_matrices(desf, hyp, 0);
    
      
  } else {

    // Set the Estimable function matrix
    gef = set_gef_matrix(des, &rank);

    // Type of sum of squares
    if (ssq == 2) {
      hyp = set_ef3_matrix(model, gef, rank, indeces, 0);
    } else {
      hyp = set_ef3_matrix(model, gef, rank, indeces, 0);
    }

    // The Least Square Means Matrix
    lsm = set_lsmeans_matrix(model, des, indeces);

    // inv(X'X) * X'
    X = set_X_matrix(des);     
       
    // C * inv(X'X) * C'
    con = set_contrast_matrices(des, hyp, 0);

  }

  // Load the variables array with the indeces of the data vectors
  if (variable_number > 0) {
    ii = 0;
    for (i=0;i<n;i++) {
      g = search_variable(xdata, variable_names[i]);
      if (g < 0) {
	fprintf(stderr, "Unable to find %s in data\n", variable_names[i]);
	exit(1);
      }
      variables[ii] = g;
      ii++;
    }
  } else {
    for (i=0;i<n;i++) {
      variables[i]=i;
    }      
  }

  // Set the result object
  res = set_result(n,
		   factor,
		   covariate,
		   model,
		   des,
		   Q,
		   reference,
		   lsmo);
  
  // We used this v vector to copy the data
  v = TYPE_ALLOC(VECTOR);
  v->data = ARRAY_ALLOC(des->n, double);

  // Iterate over the variables (rows)
  for (i=0;i<n;i++) {

    // Set the index for the variable
    g = variables[i];

    // The size of the vector defaults to the size of the design matrix
    v->n = des->n;

    // Get the variable name
    len = strlen(xdata->variables[g]);
    res->variables[i] = ARRAY_ALLOC((len+1), char);
    strncpy(res->variables[i], xdata->variables[g], len);
    res->variables[i][len] = '\0';

    if (verbose || debug) printf("%s\n", xdata->variables[g]);

    // We clean the data if there are skipped samples or NA data
    v1 = xdata->vectors[g];
    k=0;
    l=0;
    na=0;
    for (j=0;j<s;j++) {
      if (model->include[j]>0) {
	if (! isnan(v1->data[j])) {
	  v->data[k] = v1->data[j];
	  rna[k]=j-l;
	  k++; 
	} else {
	  na++;
	}
      } else {
	l++;
      }
    }

    // Resize the vector if there is NA data
    v->n -= na;

    if (na > 0) {

      // We need to create new matrices if there is NA data
      if (debug) printf("Original Matrices Modified\n");

      desna = clean_design_matrix(rna, des, na); 

      if (model->r > 0) {

	desfna = clean_design_matrix(rna, desf, na); 
	desrna = clean_design_matrix(rna, desr, na); 
	Rna = identity((model->s - model->k) - na);
	R_derivna = TYPE_ALLOC(NMATRIX);
	R_derivna->n =  1;
	R_derivna->m = ARRAY_ALLOC(R_derivna->n, MATRIX *);
	for (j = 0; j < 1; j++) {
          R_derivna->m[j] = identity((model->s - model->k) - na);
	}
	gefna = set_gef_matrix(desfna, &rank);
	if (ssq == 2) {
          hypna = set_ef3_matrix(model, gefna, rank, findeces, 0);
	} else {
          hypna = set_ef3_matrix(model, gefna, rank, findeces, 0);
	}
	lsmna = set_lsmeans_matrix(model, desna, indeces);
	Xna = set_X_matrix(desfna);
	conna = set_contrast_matrices(desfna, hypna, 0);
	mixed(v,desrna,rindeces,desfna,Rna,G,&hess,grad,G_deriv,R_derivna,desna,indeces,Xna,hypna,ndf,conna,lsmna,res,model,method,i,reference,lsmo,islog);
	free(desna);
	free(desfna);
	free(desrna);
	free(Rna);
	free(R_derivna);
	free(gefna);
	free(Xna);
	free(hypna);
	free(lsmna);
	free(conna);
	
      } else {

	gefna = set_gef_matrix(desna, &rank);
	if (ssq == 2) {
	  hypna = set_ef3_matrix(model, gefna, rank, indeces, 0);
	} else {
          hypna = set_ef3_matrix(model, gefna, rank, indeces, 0);
	}
	lsmna = set_lsmeans_matrix(model, desna, indeces);
	Xna = set_X_matrix(desna);
	conna = set_contrast_matrices(desna, hypna, 0);
	glm(v,desna,Xna,hypna,conna,lsmna,res,model,i,reference,lsmo,islog);
	free(desna);
	free(gefna);
	free(Xna);
	free(hypna);
	free(lsmna);
	free(conna);

      }

    } else {

      if (model->r > 0) {
        
	mixed(v,desr,rindeces,desf,R,G,&hess,grad,G_deriv,R_deriv,des,indeces,X,hyp,ndf,con,lsm,res,model,method,i,reference,lsmo,islog);

      } else {

	glm(v,des,X,hyp,con,lsm,res,model,i,reference,lsmo,islog);

      }

    } 

    free(v1);

  }

  mult_test_correction(res,model,Q);

  free(v);
  free(rna);
  free(rank);
  free(indeces);
  free(des);
  free(gef);
  free(X);
  free(hyp);
  free(lsm);
  free(con);

  if (model->r > 0) {
    free(findeces);
    free(rindeces);
    free(desf);
    free(desr);
    free(R);
    free(G);
    free(hess);
    free(grad);
    free(G_deriv);
    free(R_deriv);
  }

  return res;

}

int
search_variable(XDATA *xdata,
		char *variable_name)
{

  int i;
  
  for (i = 0; i < xdata->n; i++) {
    if (strcmp(variable_name, xdata->variables[i]) == 0) {
      return i;
    }
  }
  
  return (int) -1;

}

MATRIX *
set_X_matrix(MATRIX *des)
{

  int *rank, r;
  double det;
  MATRIX *Xt, *XtX, *iXtX, *X;
  
  Xt = transpose(des,0);
  XtX = multiply(Xt,des,0);
  rank = ARRAY_ALLOC(XtX->m, int);
  iXtX = g2invert(XtX,
		  XtX->m,
		  &rank,
		  &r,
		  &det,
		  1);
  X = multiply(iXtX,Xt,1);
  free(rank);
  
  return X;

}


MATRIX *
set_coef_matrix(MATRIX *XRX,
		MATRIX *XRZ,
		MATRIX *ZRX,
		MATRIX *ZRZG)
{

  MATRIX *ColU, *ColD, *Row;

  ColU = augment_cols(XRX,XRZ,1);
  ColD = augment_cols(ZRX,ZRZG,1);
  Row = augment_rows(ColU,ColD,1);

  return Row;

}

MATRIX *
set_initial_W_matrix(MATRIX *X,
		     MATRIX *Z,
		     MATRIX *Y,
		     MATRIX **bu,
		     MATRIX **varb)
{

  // This is to construct the W matrix according to Computing Gaussian
  // Likelihoods and their derivatives for general linear mixed
  // models. SIAMJ. Sci Comput.  Vol 15, No 6, page 1294-1310 November
  // 1994 This subroutine is in page 1300 section 3.2 The initial W is
  // computed with V inverse is Identity

  int *rank, r;
  double det;
  MATRIX *Xt, *XtX, *iXtX;
  MATRIX *Zt, *XCt, *XC, *b, *Xb, *R, *Rt;
  MATRIX *XCtXC, *XCtZ, *XCtR;
  MATRIX *ZtXC, *ZtZ, *ZtR;
  MATRIX *RtXC, *RtZ, *RtR;
  MATRIX *rows0, *rows1, *rows2;
  MATRIX *iW;

  // Set initial parameters for bu and varb

  Zt = transpose(Z,0);
  XCt = set_X_matrix(X);
  b = multiply(XCt,Y,0);
  Xb = multiply(X,b,0);
  R = subtract(Y,Xb,0);
  Rt = transpose(R,0);
  
  Xt = transpose(X,0);
  XtX = multiply(Xt,X,0);
  rank = ARRAY_ALLOC(XtX->m, int);
  iXtX = g2invert(XtX,
  		  XtX->m,
  		  &rank,
  		  &r,
  		  &det,
  		  1);

  *bu = fill(Z->m,1,0);
  *bu = augment_rows(b,*bu,1);
  *varb = duplicate(iXtX);

  // Set initial X*

  XC = compute_X_star_matrix(iXtX);
  XC = multiply(X,XC,0);		  
  XCt = transpose(XC,0);

  // Set initial W matrix
 
  XCtXC = multiply(XCt,XC,0);
  XCtZ = multiply(XCt,Z,0);
  XCtR =multiply(XCt,R,0);

  ZtXC = multiply(Zt,XC,0);
  ZtZ = multiply(Zt,Z,0);
  ZtR =multiply(Zt,R,0);

  RtXC = multiply(Rt,XC,0);
  RtZ = multiply(Rt,Z,0);
  RtR =multiply(Rt,R,0);

  rows0 = augment_cols(XCtXC,XCtZ,1);
  rows0 = augment_cols(rows0,XCtR,1);
  rows1 = augment_cols(ZtXC,ZtZ,1);
  rows1 = augment_cols(rows1,ZtR,1);
  rows2 = augment_cols(RtXC,RtZ,1);
  rows2 = augment_cols(rows2,RtR,1);
  
  iW = augment_rows(rows0,rows1,1);
  iW = augment_rows(iW,rows2,1);

  if (debug) {
    printf ("INITIAL W TRANSFORM MATRIX %iX%i\n\n", iW->n, iW->m);
    dump_matrix(iW);
    printf ("\n\n");
  }
  
  free(Xt);
  free(iXtX);
  free(rank);
  free(Zt);
  free(XCt);
  free(XC);
  free(Xb);
  free(R);
  free(Rt);

  return iW;
  
}

MATRIX *
compute_W_transform(MATRIX *W0a,
		    MATRIX *G,
		    MATRIX *X,
		    MATRIX *R,
		    MATRIX **bu,
		    MATRIX **varb,
		    enum METHOD method,
		    double *ml,
		    double *reml,
		    int n)
{

    // We are not using this method currently (typical
    // Newton-Raphson method used in SAS). We use instead
    //compute_EM_matrix

    // This is to construct the W0 matrix according to 
    // Computing Gaussian Likelihoods and their derivatives
    // for general linear mixed models. SIAMJ. Sci Comput.
    // Vol 15, No 6, page 1294-1310 November 1994
    // This subroutine is in page 1296-8 section 2.1 and 2.3
    // We pass the inverted R matrix to improve
    // performance
  
    int *rank, r, rr;
    double logR, l2;
    double l1, l3;
    MATRIX *cW0a, *W, *cW, *b, *u;

    // Computation of the augmented W0 sweep
    // Matrix for derivative estimation
    cW0a = duplicate(W0a);

    rank = ARRAY_ALLOC(cW0a->m, int);
    sweep1(cW0a,
	   0,
	   G->m,
	   &rank,
	   &r,
	   &l1);
    // Submatrix of augmented W0 sweep
    W = copy(cW0a, G->n, G->m, cW0a->n, cW0a->m);
    cW = duplicate(W);
    if (n == 0) {

	// Calculate components of the likelihood
	free(rank);
	rank = ARRAY_ALLOC(cW->m, int);		      
	sweep1(cW,
	       0,
	       X->m,
	       &rank,
	       &r,
	       &l3);
	       
	*bu = copy(cW,0,cW->m-1,cW->n-1,cW->m);
	b = copy(*bu,0,0,X->m,1);
	u = copy(*bu,X->m,0,(*bu)->n,1); 
	u = multiply(G,u,0);
	*bu = augment_rows(b,u,0);
	*varb = copy(cW,0,0,X->m,X->m);
	logR = trace_log(R, 0);
	l2 = cW->data[cW->n-1][cW->m-1];
	if (method == ML) {
	    *ml = l1+logR+l2+((X->n)*log(2*PI));
	} else  {
	    rr = rank_matrix(X);
	    *reml = l1+logR+l2+l3+((X->n-rr)*log(2*PI));      
	}

    }
  
    if (debug) {
	printf ("W TRANSFORM MATRIX %iX%i\n\n", W->n, W->m);
	dump_matrix(W);
	printf ("\n\n");
	printf ("SOLUTION MATRIX FOR LIKELIHOOD %iX%i\n\n", cW->n, cW->m);
	dump_matrix(cW);
	printf ("\n\n");
	if (method == ML) {
	    printf("ML = %f\n", *ml);
	} else {
	    printf("REML = %f\n", *reml);
	}
    }

    free(rank); 
    free(cW0a);
    free(cW);
    if (n == 0) {
      free(b);
      free(u);
    }
    return W;

}

MATRIX *
compute_W0_augment(MATRIX *W0,
		   MATRIX *iR,
		   MATRIX *X,
		   MATRIX *Z,
		   MATRIX *G)
{

  // This is to construct the W0 matrix according to 
  // Computing Gaussian Likelihoods and their derivatives
  // for general linear mixed models. SIAMJ. Sci Comput.
  // Vol 15, No 6, page 1294-1310 November 1994
  // This subroutine is in page 1297 section 2.3
  // We pass the inverted R matrix to improve
  // performance

  float *p;
  MATRIX *L, *I, *Lt, *ZRZ, *LtZRZ, *LtZRZL, *ILtZRZL;
  MATRIX *WZ, *LtWZ, *LtWZt;
  MATRIX *W0a;

  // Allocate memory
  p = ARRAY_ALLOC(G->n, float);

  // Calculation of Cholesky Decomposition
  L = duplicate(G);
  cholroot(L->data, L->n, p);

  // Upper Left of W0 augment
  I = identity(L->n);
  ZRZ = set_W_trans(Z,iR,Z);
  Lt = transpose(L,0);
  LtZRZ = multiply(Lt,ZRZ,0);
  LtZRZL = multiply(LtZRZ,L,1);
  ILtZRZL = add(I,LtZRZL,1);
  
  // Upper Right of W0 augment
  WZ = copy(W0, X->m, 0, X->m + Z->m, W0->m); 
  LtWZ = multiply(Lt,WZ,1);

  // Lower Left of W0 augment
  LtWZt = transpose(LtWZ,0);

  W0a = set_coef_matrix(ILtZRZL,LtWZ,LtWZt,W0);

  if (debug) {
    printf ("W0 AUGMENTED MATRIX %iX%i\n\n", W0->n, W0->m);
    dump_matrix(W0);
    printf ("\n\n");    
  }
  
  free(ZRZ);

  return W0a;

}

MATRIX *
compute_W0_transform(MATRIX *X,
		     MATRIX *Z,
		     MATRIX *iR,
		     MATRIX *Y)
{

  // This is to construct the W0 matrix according to 
  // Computing Gaussian Likelihoods and their derivatives
  // for general linear mixed models. SIAMJ. Sci Comput.
  // Vol 15, No 6, page 1294-1310 November 1994
  // This subroutine is in page 1297 section 2.3
  // We pass the inverted R matrix to improve
  // performance

  MATRIX *XRX, *XRZ, *ZRX, *ZRZ, *XRY, *ZRY;
  MATRIX *coeff, *yBord, *yBordt, *yR, *W0;
  
  // W0 Matrix (Left side)
  XRX = set_W_trans(X,iR,X);
  XRZ = set_W_trans(X,iR,Z);
  ZRX = set_W_trans(Z,iR,X);
  ZRZ = set_W_trans(Z,iR,Z);
  coeff = set_coef_matrix(XRX,XRZ,ZRX,ZRZ);

  // W0 Matrix (Right side)
  XRY = set_W_trans(X,iR,Y);
  ZRY = set_W_trans(Z,iR,Y);
  yBord = set_y_border(XRY,ZRY);  
  yBordt = transpose(yBord,0);
  yR = set_W_trans(Y,iR,Y);
  
  // W0 Matrix
  W0 = set_coef_matrix(coeff,yBord,yBordt,yR);

  if (debug) {
    printf ("W0 MATRIX %iX%i\n\n", W0->n, W0->m);
    dump_matrix(W0);
    printf ("\n\n");    
  }

  return W0;

}

MATRIX *
set_W_trans(MATRIX *X,
	    MATRIX *iR,
	    MATRIX *Y)
{

  MATRIX *Xt, *XtiR, *XtiRY;
  
  Xt = transpose(X,0);
  XtiR = multiply(Xt,iR,0);
  XtiRY = multiply(XtiR,Y,0);
  free(Xt);
  free(XtiR);
  
  return XtiRY;

}

MATRIX *
set_y_border(MATRIX *XRY,
		MATRIX *ZRY)
{

  MATRIX *X;

  X = augment_rows(XRY,ZRY,1);

  return X;

}

MATRIX *
compute_subvariance_matrix(MATRIX *iR,
			   MATRIX *Z,
			   MATRIX *G)
{

  // This is to construct the W0 matrix according to 
  // Computing Gaussian Likelihoods and their derivatives
  // for general linear mixed models. SIAMJ. Sci Comput.
  // Vol 15, No 6, page 1294-1310 November 1994
  // This subroutine is in page 1301 section 3.3
  // We pass the inverted R matrix to improve
  // performance  


  int *rank, r;
  double det;
  MATRIX *iG, *ZRZ, *iM, *M;
  
  rank = ARRAY_ALLOC(G->m, int);
  iG = g2invert(G,
		G->m,
		&rank,
		&r,
		&det,
		0);
  ZRZ = set_W_trans(Z,iR,Z);
  iM = add(iG,ZRZ,0);
  M = g2invert(iM,
	       iM->m,
	       &rank,
	       &r,
	       &det,
	       1);

  if (debug) {
    printf ("M MATRIX %iX%i\n\n", M->n, M->m);
    dump_matrix(M);
    printf ("\n\n");    
  }

  free(rank);
  free(iG);
  free(ZRZ);

  return M;

}

MATRIX *
compute_variance_matrix(MATRIX *iR,
			MATRIX *Z,
			MATRIX *M)
{

  // This is to construct the W0 matrix according to 
  // Computing Gaussian Likelihoods and their derivatives
  // for general linear mixed models. SIAMJ. Sci Comput.
  // Vol 15, No 6, page 1294-1310 November 1994
  // This subroutine is in page 1301 section 3.3
  // We pass the inverted R matrix to improve
  // performance  

  MATRIX *Zt, *ZMZ, *ZMZiR;
  MATRIX *I, *IZMZiR, *V;
  
  Zt = transpose(Z,0);
  ZMZ = set_W_trans(Zt,M,Zt);
  ZMZiR = multiply(ZMZ,iR,0);
  I = identity(ZMZ->n);
  IZMZiR = subtract(I,ZMZiR,1);
  V = multiply(iR,IZMZiR,0);

  free(Zt);
  free(ZMZ);
  free(IZMZiR);

  if (debug) {
    printf ("V MATRIX %iX%i\n\n", V->n, V->m);
    dump_matrix(V);
    printf ("\n\n");    
  }

  return V;

}

MATRIX *
compute_X_star_matrix(MATRIX *V)

{

  // This is to construct the W0 matrix according to 
  // Computing Gaussian Likelihoods and their derivatives
  // for general linear mixed models. SIAMJ. Sci Comput.
  // Vol 15, No 6, page 1294-1310 November 1994
  // This subroutine is in page 1299 section 3.1

  float *p;
  MATRIX *C;

  // Allocate memory
  p = ARRAY_ALLOC(V->n, float);

  C = duplicate(V);
  cholroot(C->data, C->n, p);

  return C;

}

void
compute_derivatives(MATRIX *X,
		    MATRIX *Z,
		    MATRIX *iR,
		    MATRIX *G,
		    MATRIX *Y,
		    MATRIX **hess,
		    MATRIX **grad,
		    MATRIX *W,
		    MATRIX *bu,
		    MATRIX *varb,
		    NMATRIX *G_deriv,
		    NMATRIX *R_deriv,
		    enum METHOD method,
		    int n)
{
  // This is to compute derivaties for Newton-Raphson algorithm 
  // to solve mixed model equations according to
  // Computing Gaussian Likelihoods and their derivatives
  // for general linear mixed models. SIAMJ. Sci Comput.
  // Vol 15, No 6, page 1294-1310 November 1994
  // This subroutine is in page 1300 section 3.2

  int i, j;
  double g1, g2, g3, h1, h2, h3, sum;
  MATRIX *WXX, *WZX, *WrX;
  MATRIX *WXZ, *WZZ, *WrZ;
  MATRIX *WXr, *WZr, *Wrr;  
  MATRIX *H1rs, *H2r, *H2rt, *H2s, *H2rs, *H3r, *H3s, *H3rs;
  MATRIX *g1r, *g2r;
  MATRIX *H2rs2, *H3rs2;
  MATRIX *M, *V, *XC, *XCt, *b, *r, *rt;
  MATRIX *ZGZ, *Zt;

  if (n > 1) {
    sum = 0.0;
    for (i=0;i<G->n;i++) {
      sum += fabs(G->data[i][i]);
    } 
    if (sum < EPS) {
      M = fill(G->n,G->n,0);
    } else {
      M = compute_subvariance_matrix(iR, Z, G);
    }
    V = compute_variance_matrix(iR, Z, M);
 
  } else {
    V = identity(X->n);
  }

  // Compute X*

  XC = compute_X_star_matrix(varb);
  XC = multiply(X, XC, 0);
  XCt = transpose(XC,0);

  // Compute r = (Y - Xb)

  Zt = transpose(Z,0);
  b = copy(bu,0,0,X->m,bu->m);
  b = multiply(X,b,0);
  r = subtract(Y,b,0);
  rt = transpose(r,0);  

  // Compute derivatives for parameters in G using the W matrix (Page
  // 1300)
 
  WXX = copy(W,0,0,X->m,X->m);
  WZX = copy(W,X->m,0,X->m+Z->m,X->m);
  WrX = copy(W,X->m+Z->m,0,W->n,X->m);

  WXZ = copy(W,0,X->m,X->m,X->m+Z->m);
  WZZ = copy(W,X->m,X->m,X->m+Z->m,X->m+Z->m);
  WrZ = copy(W,X->m+Z->m,X->m,W->n,X->m+Z->m);

  WXr = copy(W,0,X->m+Z->m,X->m,W->m);
  WZr = copy(W,X->m,X->m+Z->m,X->m+Z->m,W->m); 
  Wrr = copy(W,X->m+Z->m,X->m+Z->m,W->n,W->m);
  
  for (i = 0; i < (*hess)->m-1; i++) {

    H2r = multiply(WXZ,G_deriv->m[i],0);
    H2r = multiply(H2r,WZr,0);
    H2rt = transpose(H2r,0);
    H3r = multiply(WXZ,G_deriv->m[i],0);
    H3r = multiply(H3r,WZX,0);
    g1r = multiply(WZZ,G_deriv->m[i],0);    
    g2r = multiply(WrZ,G_deriv->m[i],0);
    g2r = multiply(g2r,WZr,0);    
    g1 = trace(g1r,0);
    g2 = (-1) * g2r->data[0][0];
     
    for (j = 0; j < (*hess)->m-1; j++) {
  
      H2s = multiply(WXZ,G_deriv->m[j],0);
      H2s = multiply(H2s,WZr,0);
      H3s = multiply(WXZ,G_deriv->m[j],0);
      H3s = multiply(H3s,WZX,0);
      H1rs = multiply(g1r, WZZ, 0);
      H1rs = multiply(H1rs, G_deriv->m[j],0);      
      H2rs = multiply(WrZ,G_deriv->m[i],0); 
      H2rs = multiply(H2rs,WZZ,0);
      H2rs = multiply(H2rs,G_deriv->m[j],0);
      H2rs = multiply(H2rs,WZr,0);
      H2rs = multiply_number(H2rs,2,0);
      H3rs = multiply(WXZ,G_deriv->m[i],0);
      H3rs = multiply(H3rs,WZZ,0);
      H3rs = multiply(H3rs,G_deriv->m[j],0);
      H3rs = multiply(H3rs,WZX,0);
      H3rs = multiply_number(H3rs,2,0);
      H2rs2 = multiply(H2rt,H2s,0);
      H2rs2 = multiply_number(H2rs2,2,0);
      H2rs2 = subtract(H2rs,H2rs2,0);
      H3rs2 = multiply(H3r,H3s,0);
      H3rs2 = subtract(H3rs,H3rs2,0);
      
      g3 = (-1) * trace(H3r,0);
      h1 = (-1) * trace(H1rs,0);
      h2 = H2rs2->data[0][0];
      h3 = trace(H3rs2,0);
            
      if (method == ML && n>1) {
        (*grad)->data[i][0] = g1 + g2;
        (*hess)->data[i][j] = h1 + h2; 
      } else if (method == REML && n>1) {
        (*grad)->data[i][0] = g1 + g2 + g3;
        (*hess)->data[i][j] = h1 + h2 + h3; 
      } else {
        (*grad)->data[i][0] = g2;
        (*hess)->data[i][j] = h1 + h3;      
      }

    }
    
  }
  
  // Compute derivatives for parameters in R using M and V (computed
  // above) (Page 1301)
  
  for (i = (*hess)->m-1; i < (*hess)->m; i++) {
    
    H2r = multiply(XCt,V,0);
    H2r = multiply(H2r,V,0);
    H3r = multiply(H2r,XC,0);
    H2r = multiply(H2r,r,0);
    H2rt = transpose(H2r,0);
    g2r = multiply(rt,V,0);
    g2r = multiply(g2r,V,0);
    g2r = multiply(g2r,r,0);        
    g1 = trace(V,0);
    g2 = (-1) * g2r->data[0][0]; 
 
    for (j = (*hess)->m-1; j < (*hess)->m; j++) {
        
      H2s = H2r;
      H3s = H3r;
      H1rs = multiply(V, V, 0);
      H2rs = multiply(rt,V,0);
      H2rs = multiply(H2rs,V,0);
      H2rs = multiply(H2rs,V,0);
      H2rs = multiply(H2rs,r,0);
      H2rs = multiply_number(H2rs,2,0);
      H3rs = multiply(XCt,V,0);
      H3rs = multiply(H3rs,V,0);
      H3rs = multiply(H3rs,V,0);
      H3rs = multiply(H3rs,XC,0);
      H3rs = multiply_number(H3rs,2,0);
      H2rs2 = multiply(H2rt,H2s,0);
      H2rs2 = multiply_number(H2rs2,2,0);
      H2rs2 = subtract(H2rs,H2rs2,0);
      H3rs2 = multiply(H3r,H3s,0);
      H3rs2 = subtract(H3rs,H3rs2,0);
      
      g3 = (-1) * trace(H3r,0);
      h1 = (-1) * trace(H1rs,0);
      h2 = H2rs2->data[0][0];
      h3 = trace(H3rs2,0);
            
      if (method == ML && n>1) {
        (*grad)->data[i][0] = g1 + g2;
        (*hess)->data[i][j] = h1 + h2; 
      } else if (method == REML && n>1) {
        (*grad)->data[i][0] = g1 + g2 + g3;
        (*hess)->data[i][j] = h1 + h2 + h3; 
      } else {
        (*grad)->data[i][0] = g2;
        (*hess)->data[i][j] = h1 + h3;      
      } 

    }
  
  }
  
  // Compute cross derivatives for parameters in R and G using M and V
  // (computed above) (Page 1302-3)
  
  for (i = 0; i < (*hess)->m-1; i++) {

    for (j = (*hess)->m-1; j < (*hess)->m; j++) {
    
      ZGZ = multiply(Z,G_deriv->m[i],0);
      ZGZ = multiply(ZGZ,Zt,0);
      
      H2r = multiply(XCt,V,0);
      H2s = multiply(H2r,V,0);
      H2s = multiply(H2s,r,0);
      H3s = multiply(H2r,V,0);
      H3s = multiply(H3s,XC,0);
      H2r = multiply(H2r,ZGZ,0);
      H2r = multiply(H2r,V,0);
      H3r = multiply(H2r,XC,0);
      H2r = multiply(H2r,r,0);
      H2rt = transpose(H2r,0);
      H1rs = multiply(V, ZGZ, 0);
      H1rs = multiply(H1rs, V, 0);
      H2rs = multiply(rt,V,0);
      H2rs = multiply(H2rs,ZGZ,0);
      H2rs = multiply(H2rs,V,0);
      H2rs = multiply(H2rs,V,0);
      H2rs = multiply(H2rs,r,0);
      H2rs = multiply_number(H2rs,2,0);
      H3rs = multiply(XCt,V,0);
      H3rs = multiply(H3rs,ZGZ,0);
      H3rs = multiply(H3rs,V,0);
      H3rs = multiply(H3rs,V,0);
      H3rs = multiply(H3rs,XC,0);
      H3rs = multiply_number(H3rs,2,0);
      H2rs2 = multiply(H2rt,H2s,0);
      H2rs2 = multiply_number(H2rs2,2,0);
      H2rs2 = subtract(H2rs,H2rs2,0);
      H3rs2 = multiply(H3r,H3s,0);
      H3rs2 = subtract(H3rs,H3rs2,0);

      h1 = (-1) * trace(H1rs,0);
      h2 = H2rs2->data[0][0];
      h3 = trace(H3rs2,0);
            
      if (method == ML && n>1) {
        (*hess)->data[i][j] = h1 + h2; 
	(*hess)->data[j][i] = h1 + h2; 
      } else if (method == REML && n>1) {
        (*hess)->data[i][j] = h1 + h2 + h3;
	(*hess)->data[j][i] = h1 + h2 + h3;
      } else {
        (*hess)->data[i][j] = h1 + h3;   
	(*hess)->data[j][i] = h1 + h3;   
      } 

    }
  
  }
   
  if (debug) {
    printf ("HESSIAN MATRIX %iX%i\n\n", (*hess)->n, (*hess)->m);
    dump_matrix(*hess);
    printf ("\n\n");   
    printf ("GRADIENT MATRIX %iX%i\n\n", (*grad)->n, (*grad)->m);
    dump_matrix(*grad);
    printf ("\n\n");   
  }

}

MATRIX *
get_Newton_direction_matrix(MATRIX **hess,
			    MATRIX *grad,
			    int p)
{

  int *rank, r;
  double det;
  MATRIX *ihess, *N;
  

  // Modification to Newton-Raphson done in
  // SAS ??? were commented out

  //det = 0.0;
  //for (i=0;i<(*hess)->n;i++) {
  //  det +=fabs((*hess)->data[i][i]);
  //}
  //det = det/(*hess)->n;
  
  //r = rank_matrix(*hess);
  //if (r<(*hess)->n+1 || p>0) {
  //  printf("Hessian is singular\n");
  //  ridge = identity((*hess)->n);
  //  ridge = multiply_number(ridge,2,0);
  //  *hess = add(*hess,ridge,0);
  //  dump_matrix(*hess);
  //}

  rank = ARRAY_ALLOC((*hess)->m, int);
  ihess = g2invert(*hess,
		   (*hess)->m,
		   &rank,
		   &r,
		   &det,
		   0);
  N = multiply(ihess,grad,0);

 // for (i=0; i< N->n; i++) {  
 //   if (N->data[i][0]<0) {
 //      N->data[i][0] = 0;
 //   }
 // }
  
  free(rank);
  free(ihess);

  return N;

}

void
Newton_Raphson(MATRIX *X,
	       MATRIX *Z,
	       MATRIX *iR,
	       MATRIX *G,
	       MATRIX *R,
	       int **rindeces,
	       MATRIX *Y,
	       MATRIX *parameters,
	       MATRIX **hess,
	       MATRIX **grad,
	       MATRIX **bu,
	       MATRIX **varb,
	       NMATRIX *G_deriv,
	       NMATRIX *R_deriv,
	       enum METHOD method,
	       double *ml,
	       double *reml)
{

  int i, j, c=0, *rank, r, maxiter = 10, ridge = 0;
  double det, sum, oml, oreml;
  MATRIX *XC, *b, *Xb, *YXb, *W0, *Wa, *W;
  MATRIX *iparameters;

  iparameters = duplicate(parameters);
  oml = *ml;
  oreml = *reml;
  // Loop for Newton - Raphson

  for (i = 2; i < maxiter; i++) {

    // Compute X star  

    XC = compute_X_star_matrix(*varb);
    XC = multiply(X,XC,0);
  
    // Compute W for derivaties
    b = copy(*bu,0,0,X->m,(*bu)->m);
    Xb = multiply(X,b,0);   
    YXb = subtract(Y,Xb,0);
    W0 = compute_W0_transform(XC,Z,iR,YXb);
    Wa = compute_W0_augment(W0,iR,X,Z,G);
    W = compute_W_transform(Wa,G,X,R,&*bu,&*varb,method,&*ml,&*reml,1);
    compute_derivatives(X,Z,iR,G,Y,&*hess,&*grad,W,*bu,*varb,G_deriv,R_deriv,method,i); 

    // Get new parameters    
    parameters = get_Newton_direction_matrix(&*hess, *grad, 0);
    parameters = subtract(iparameters, parameters, 0);

    while (ridge == 0 && c<100) {
  
	for (j=0; j< parameters->n; j++) {  
	    if (parameters->data[j][0]<0) {
		parameters->data[j][0] = 0;
	    }
	}
    
	// Update G and R matrices
	update_G_matrix(&G,parameters,rindeces);
	update_R_matrix(&R,parameters);  

	// We invert the r Matrix only
	// once to improve performance
	rank = ARRAY_ALLOC(R->m, int);
	iR = g2invert(R,
		      R->m,
		      &rank,
		      &r,
		      &det,
		      0);

	// Likelihood Calculations
	W0 = compute_W0_transform(X,Z,iR,Y);
	Wa = compute_W0_augment(W0,iR,X,Z,G);
	W = compute_W_transform(Wa,G,X,R,&*bu,&*varb,method,&*ml,&*reml,0);
	if (*reml < oreml) {
	    ridge = 1;
	} else {
	    c++;
	    iparameters = duplicate(parameters);
	    oreml = *reml;
	    parameters = get_Newton_direction_matrix(&*hess, *grad, 1);
	    parameters = subtract(iparameters, parameters, 1);
	    ridge = 1;
	}

    }   
    // Test if we converged
    sum = 0;
    for (j = 0; j < parameters->n; j++) {
      sum += pow(parameters->data[j][0] - iparameters->data[j][0], 2);
    }
    if (sum < EPS) {
         printf("SUM is %f\n",sum);
         return;
    }
	
    // Copy the new parameters
    iparameters = duplicate(parameters);

  }   

}

MATRIX * 
compute_EM_matrix(MATRIX *X,
	          MATRIX *Z,
	          MATRIX *G,
	          MATRIX *Y,
		  MATRIX *parameters) 
{

  // Method for minimizing REML and ML from MATLAB 
  // MATLAB algorithm mixed.m for solving Henderson's mixed model equations
  // Tecnhical Report, Institute of Measurement Science, Slovak Academy of Sciences, Bratislava, Deecember 2001, page 4
  
  int *rank,r;
  double det;
  MATRIX *Xt, *XtX, *XtZ, *Zt, *ZtX, *ZtZ, *id, *V, *XtY, *ZtY, *XZY, *C, *iC;
  MATRIX *bu;

  Xt = transpose(X,0);
  XtX = multiply(Xt,X,0);
  XtZ = multiply(Xt,Z,0);
  XtZ = multiply(XtZ,G,0);
  Zt = transpose(Z,0);
  ZtX = multiply(Zt,X,0);
   
  id = identity(G->n);
 
  id = multiply_number(id, parameters->data[parameters->n-1][0],0);
  ZtZ = multiply(Zt,Z,0);
  ZtZ = multiply(ZtZ,G,0);
  V = add(id, ZtZ,1);

  C = set_coef_matrix(XtX,XtZ,ZtX,V);
  XtY = multiply(Xt,Y,0);
  ZtY = multiply(Zt,Y,0);
  XZY = augment_rows(XtY,ZtY,1);
  
  rank = ARRAY_ALLOC(C->m, int);
  iC = g2invert(C,
		C->m,
		&rank,
		&r,
		&det,
		0);
  bu = multiply(iC,XZY,1);
  
  free(rank);
  free(Xt);
  free(Zt);
  free(C);

  return bu;
  	  
}

void
Expected_Maxim(MATRIX *X,
	       MATRIX *Z,
	       MATRIX *iR,
	       MATRIX *G,
	       MATRIX *R,
	       int **rindeces,
	       MATRIX *Y,
	       MATRIX *parameters,
	       MATRIX **hess,
	       MATRIX **grad,
	       MATRIX **bu,
	       MATRIX **varb,
	       NMATRIX *G_deriv,
	       NMATRIX *R_deriv,
	       enum METHOD method,
	       double *ml,
	       double *reml)
{


  // Method for minimizing REML and ML from MATLAB 
  // MATLAB algorithm mixed.m for solving Henderson's mixed model equations
  // Tecnhical Report, Institute of Measurement Science, Slovak Academy of Sciences, Bratislava, Deecember 2001
  
  int i, j, *rank, r, rx, maxiter = 400, d0=0;
  double det, sum, k, lr, oreml, oml;
  MATRIX *b, *u, *cu, *cut, *W, *W0, *Wa, *iW, *ciW, *id, *Zt, *ZtZG, *M;
  MATRIX *Xb, *Zu, *Ys, *Yt;
  MATRIX *iparameters, *tbu;
  
  oreml = *reml;
  oml = *ml;
  iparameters = duplicate(parameters);
   
  // Loop for EM

  for (i = 1; i < maxiter; i++) {
  
    // Compute solutions for derivaties page 4 in MATLAB algorithm mixed.m for solving Henderson's mixed model equations    
   
    tbu = compute_EM_matrix(X,Z,G,Y,parameters);
    b = copy(tbu,0,0,X->m,tbu->m);  
    Xb = multiply(X,b,0);      
    u = copy(tbu,X->m,0,tbu->n,tbu->m);
    u = multiply(G,u,0);
    Zu = multiply(Z,u,0);
    Ys = add(Xb,Zu,0);
    Ys = subtract(Y,Ys,0);
    Yt = transpose(Y,0);
    Yt = multiply(Yt,Ys,0);
    id = identity(G->n);
    Zt = transpose(Z,0);
    ZtZG = multiply(Zt,Z,0);
    ZtZG = multiply(ZtZG,G,0);
    if (method == REML) {
      id = identity(X->n);
      M = set_X_matrix(X);
      M = multiply(X,M,0);
      M = subtract(id,M,0);
      ZtZG = set_W_trans(Z,M,Z);
      ZtZG = multiply(ZtZG,G,0);
    }
   
    id = identity(G->n);
    W = multiply_number(id,parameters->data[parameters->n-1][0],0);
    W = add(W,ZtZG,0);
    rank = ARRAY_ALLOC(W->m, int);
    iW = g2invert(W,
	  	  W->m,
		  &rank,
		  &r,
		  &det,
		  0);
    iW = multiply_number(iW,parameters->data[parameters->n-1][0],0);
    ciW = duplicate(iW);

    free(rank);
    rank = ARRAY_ALLOC(ciW->m, int);
    sweep(ciW,
	  0,
	  ciW->m,
	  &rank,
	  &r,
	  &lr);

    for (j=0;j<parameters->n-1;j++) {
      ciW = copy(iW,rindeces[j][0],rindeces[j][0],rindeces[j][1],rindeces[j][1]);
      k = rindeces[j][1]-rindeces[j][0]-trace(ciW,0);
      cu = copy(u,rindeces[j][0],0,rindeces[j][1],u->m);
      cut = transpose(cu,0);
      cut = multiply(cut,cu,0);
      if (k > 100*DBL_EPSILON) {
        parameters->data[j][0] = cut->data[0][0]/k;
      } else {
        parameters->data[j][0] = min_matrix(parameters);
	d0++;
      }
    }
    
    if (method == REML) {
      rx = rank_matrix(X);
      parameters->data[parameters->n-1][0] = Yt->data[0][0]/((X->n)-rx);
    } else {
      parameters->data[parameters->n-1][0] = Yt->data[0][0]/(X->n);
    }   
    
    // Update G and R matrices
    update_G_matrix(&G,parameters,rindeces);
    update_R_matrix(&R,parameters);  
     
    // We invert the r Matrix only
    // once to improve performance
    rank = ARRAY_ALLOC(R->m, int);
    iR = g2invert(R,
		  R->m,
		  &rank,
		  &r,
		  &det,
		  0);

    // Likelihood Calculations
   
    if (method == REML) {
      *reml = ((X->n - rx)*log(2*PI*parameters->data[parameters->n-1][0])-lr+(X->n-rx));
    } else {
      *ml = ((X->n)*log(2*PI*parameters->data[parameters->n-1][0])-lr+X->n);
    }
       
    // Test if we converged
    sum = 0.0;
    for (j = 0; j < parameters->n; j++) {
      sum += pow(fabs(parameters->data[j][0] - iparameters->data[j][0]), 2);
    }

    if (sum < TOL) {
      
      // Final REML, ML calculations and parameters estimations
      W0 = compute_W0_transform(X,Z,iR,Y);
      Wa = compute_W0_augment(W0,iR,X,Z,G);
      W = compute_W_transform(Wa,G,X,R,&*bu,&*varb,method,&*ml,&*reml,0);
      return;

    }
	
    W0 = compute_W0_transform(X,Z,iR,Y);
    Wa = compute_W0_augment(W0,iR,X,Z,G);
    W = compute_W_transform(Wa,G,X,R,&*bu,&*varb,method,&*ml,&*reml,0);
    // Stop if REML increases;
    if (method == REML) {
      if ((*reml) > oreml) {
        update_G_matrix(&G,iparameters,rindeces);
        update_R_matrix(&R,iparameters);  
        i = maxiter; 
      } else {
        oreml = (*reml);
      }
    } else { 
      if ((*ml) > oml) {
        update_G_matrix(&G,iparameters,rindeces);
        update_R_matrix(&R,iparameters);  
        i = maxiter; 
      } else {
        oml = (*ml);
      }
    }

    // Copy the new parameters
    iparameters = duplicate(parameters);

  }  
  
  // Final REML, ML calculations and parameters estimations
  W0 = compute_W0_transform(X,Z,iR,Y);
  Wa = compute_W0_augment(W0,iR,X,Z,G);
  W = compute_W_transform(Wa,G,X,R,&*bu,&*varb,method,&*ml,&*reml,0);
  
  free(rank);
  free(b);
  free(u);
  free(cu);
  free(cut);
  free(W);
  free(iW);
  free(ciW);
  free(id);
  free(Zt);
  free(ZtZG);
  if (method == REML) {
    free(M);
  }
  free(Xb);
  free(Zu);
  free(Ys);
  free(Yt);
  
}

NMATRIX *
set_contrast_matrices(MATRIX *des,
		      NMATRIX *hyp,
                      int type)
{

  int i;
  int *rank, r;
  double det;
  MATRIX *Xt, *XtX, *iXtX;
  MATRIX *Ct, *CiXtX, *CiXtXCt;
  NMATRIX *con;  

  // Memory allocation
  con = TYPE_ALLOC(NMATRIX);

  con->n = hyp->n;
  con->m = ARRAY_ALLOC(con->n, MATRIX *);

  // The design matrix
  Xt = transpose(des,0);
  XtX = multiply(Xt,des,0);
  rank = ARRAY_ALLOC(XtX->m, int);
  iXtX = g2invert(XtX,
		  XtX->m,
		  &rank,
		  &r,
		  &det,
		  1);
  free(Xt);
  free(rank);

  // The contrast matrices
  for (i=0;i<con->n;i++) {
    Ct = transpose(hyp->m[i],0);
    CiXtX = multiply(hyp->m[i],iXtX,0);
    CiXtXCt = multiply(CiXtX, Ct, 1);
    if (type == 0) {
      rank = ARRAY_ALLOC(CiXtXCt->m, int);
      con->m[i] = g2invert(CiXtXCt,
    	  		   CiXtXCt->m,
    			   &rank,
    			   &r,
    			   &det,
    			   1);
      free(rank);
    } else {
      con->m[i] = CiXtXCt;
    }
  }

  free(iXtX);

  return con;

}


MATRIX *
sas_ortogonalize(MATRIX *ef,
		 MATRIX *in)
{

  int i, j, ii;
  double *w, **v, ab, bb, r;
  VECTOR *wv;
  MATRIX *m, *wm;

  if (ef->n != in->n) {
    exit_on_error("not the same number of columns.");
  }

  // Memory Allocation and initialization
  // of the orthonormal matrix
  m = TYPE_ALLOC(MATRIX); 
  m->n = ef->n;
  m->m = ef->m;
  m->data = ARRAY_ALLOC(ef->n, double *);
  for (i=0;i<ef->n;i++) {
    m->data[i] = ARRAY_ALLOC(ef->m, double);
    for (j=0;j<ef->m;j++) {
      m->data[i][j] = ef->data[i][j];
    }
  }

  // Memory allocation for the SVD method
  // We always assume that m is a subset of n 
  w = ARRAY_ALLOC(in->n, double);
  v = ARRAY_ALLOC(in->m, double *);
  for (i=0;i<in->m;i++) {
    v[i] = ARRAY_ALLOC(in->m, double);
  }

  svdcmp(in->data, in->n, in->m, w, v);

  // Convert arrays to objects
  wv = TYPE_ALLOC(VECTOR); 
  wv->n = in->m;
  wv->data = ARRAY_ALLOC(in->m, double);
  for (i=0;i<in->m;i++) {
    wv->data[i] = w[i];
  }

  wm = TYPE_ALLOC(MATRIX); 
  wm = diag(wv);
  in = multiply(in,wm,1);

  for (i=0;i<ef->m;i++) {
    for (ii=0;ii<in->m;ii++) {
      ab = dot(ef,in,i,ii);
      bb = dot(in,in,ii,ii);
      if (bb!=0) {
	r = ab/bb;
      } else {
	r = 0;
      }
      for (j=0;j<ef->n;j++) {
	m->data[j][i] -= (r * in->data[j][ii]);
	if ((double)(fabs(m->data[j][i]))<EPS) {
	  m->data[j][i] = 0;
	}
      }
    }
  }

  free(v);
  free(w);
  free(in);

  return m;

}

NMATRIX *
set_ef3_matrix(MODEL *model,
	       MATRIX *gef,
               int *rank,
               int **indeces,
	       int skip)
{
	       
  int i, ii, j, k, l, m, n, o, d, dof;
  MATRIX *ef, *in;
  NMATRIX *ef3;

  // Memory allocation
  ef3 = TYPE_ALLOC(NMATRIX);

  // If there are random effects we need to find the length of the
  // matrix otherwise this length will be the number of terms plus 1
  // (the complete model)
  n = 0;
  for (i=0;i<model->t;i++) {
    if (model->randdep[i] == 0 || skip > 0) {
      n++;
    }
  }
  ef3->n = n + 1;
  ef3->m = ARRAY_ALLOC(ef3->n, MATRIX *);
  // Set the matrices for the terms
  ii = 0;
  for (i=0;i<model->t;i++) {
    if (model->randdep[i] == 0 || skip > 0) {
      dof = 0;
      n = gef->m;
      m = indeces[ii][1] - indeces[ii][0];
      ef = mem_allocate_matrix(n, m);
      for (j=indeces[ii][0];j<indeces[ii][1];j++) {
	if (rank[j]>0) {
	  for (k=0;k<ef->n;k++) {
	    ef->data[k][dof] = gef->data[k][j];
	  }
	  dof++;
	}
      }
      ef->m = dof;
      dof = 0;
      if (model->numdep[i]>0) {
	m = gef->n - 1;
	in = mem_allocate_matrix(n, m);
	for (j=0;j<model->numdep[i];j++) {
	  d = model->dep[i][j];
	  if (model->randdep[d] == 0 || skip > 0) {
	    o = 0;
	    for (k=0;k<d;k++) {
	      o += model->randdep[k];
	    }
	    for (k=indeces[d-o][0];k<indeces[d-o][1];k++) {
	      if (rank[k]>0) {
		for (l=0;l<in->n;l++) {
		  in->data[l][dof] = gef->data[l][k];
		}
		dof++;
	      }
	    }
	  }
	}
	printf("HELLO\n");
	dump_matrix(in);
	if (dof > 0) {
	  in->m = dof;
	  ef = sas_ortogonalize(ef,in);
	}
      }
      // Copy the data to the object
      ef3->m[ii] = mem_allocate_matrix(ef->m, ef->n);
      for(j=0;j<ef->m;j++) {
	for(k=0;k<ef->n;k++) {
	  ef3->m[ii]->data[j][k] = ef->data[k][j];
	}
      }
      ii++;
    }
  }

  // Set the matrix for the model
  ef3->m[ef3->n-1] = mem_allocate_matrix(ef->n-1, ef->n);
  for(i=0;i<ef->n-1;i++) {
    for(j=0;j<ef->n;j++) {
      if (i==j-1) {
	ef3->m[ef3->n-1]->data[i][j] = 1;
      } else {
	ef3->m[ef3->n-1]->data[i][j] = 0;
      }
    }
  }
  free(ef);

  if (debug) {
    printf ("TYPE III ESTIMABLE FUNCTIONS MATRICES\n\n");
    for (i=0;i<ef3->n;i++) {
      if (i == ef3->n-1) {
	printf ("MODEL MATRIX %iX%i\n\n", ef3->m[i]->n, ef3->m[i]->m);
      } else {
	printf ("TERM %i MATRIX %iX%i\n\n", i, ef3->m[i]->n, ef3->m[i]->m);
      }
      dump_matrix(ef3->m[i]);
      printf ("\n");    
    }
    printf("\n");
  }

  return ef3;
  
}

MATRIX *
set_gef_matrix(MATRIX *x,
	       int **rank)
{

  int r, i;
  double det;
  MATRIX *m, *xt, *xtx, *ixtx;

  // General Form of Estimable
  // Functions according to SAS

  xt = transpose(x,0);
  xtx = multiply(xt,x,0);
  *rank = ARRAY_ALLOC(xtx->m, int);
  ixtx = g2invert(xtx,
		  xtx->m,
		  rank,
		  &r,
		  &det,
		  0);
  m = multiply(ixtx,xtx,0);
  m = transpose(m,1);

  if (debug) {
    printf ("GEN. ESTIMABLE FUNCTIONS (GEF) MATRIX %iX%i\n\n", m->n, m->m);
    dump_matrix(m);
    printf ("\n\n");    
    printf("RANK FOR GEF MATRIX:\n[ ");
    for (i=0;i<xtx->m;i++) {
      printf("%10i", (*rank)[i]);
    }
    printf(" ]\n\n");
  }

  free(xt);
  free(xtx);
  free(ixtx);
  
  return m;

}

MATRIX *
set_random_matrix(FACTOR *factor,
		  COVARIATE *covariate,
		  MODEL *model,
		  int ***rindeces)
{

  int i, j, c, idx;
  int id, ty;
  MATRIX *m, *m1, *m2;

  int first = 0;
  int cc = 0;

  for (i=0;i<model->i;i++) {
    if (model->randinter[i] > 0) {
      cc++;
    }
  }
  // Allocate memory
  *rindeces = ARRAY_ALLOC(model->r+cc, int *);
  model->rzeros = ARRAY_ALLOC(model->r+cc, int *);
  model->rvalid = ARRAY_ALLOC(model->r+cc, int *);
  for (i=0;i<model->r+cc;i++) {
    (*rindeces)[i] = ARRAY_ALLOC(2, int); 
  }

  idx = 0;

  // We iterate over the random factors
  for (i=0;i<model->r;i++) {
    (*rindeces)[i][0] = idx;
    id = model->rand[i]; 
    ty = model->tyfixed[i];
    m1 = design(factor, covariate, model, id, ty);
    if (first) {
      m = augment_cols(m, m1, 1);
    } else {
      m = duplicate(m1);
      first++;
    }
    idx = m->m;
    (*rindeces)[i][1] = idx;
    if (ty>0) {
      c = 1; 
    } else {
      c = factor->classes[id];
    }
    model->rzeros[i] = ARRAY_ALLOC(c, int);
    model->rvalid[i] = ARRAY_ALLOC(c, int);
    for (j=0;j<c;j++) {
      model->rzeros[i][j] = j;
      model->rvalid[i][j] = 1;
    }
  }  

  cc =  model->r;

  // We continue with the interactions
  for (i=0;i<model->i;i++) {
    if (model->randinter[i] > 0) {
      m1 = fill(model->s-model->k, 1, 1);
      for (j=0;j<model->numinter[i];j++) {
	id = model->inter[i][j];
	ty = model->tyinter[i][j];
	m2 = design(factor, covariate, model, id, ty);
	m1 = cross(m1, m2, 1);
      }
      c = m1->m;
      model->rzeros[cc] = ARRAY_ALLOC(c, int);
      model->rvalid[cc] = ARRAY_ALLOC(c, int);
      m1 = clean_cross(m1, &model->rzeros[cc]);
      for (j=0;j<c;j++) {
	if (model->rzeros[cc][j] < 0) {
	    model->rvalid[cc][j] = 0;
	} else {
	    model->rvalid[cc][j] = 1;	  
	}
      }
      (*rindeces)[cc][0] = idx;
      m = augment_cols(m, m1, 1);
      idx = m->m;
      (*rindeces)[cc][1] = idx;
      cc++;
    }
  }

  if (debug) {
    c = 0;
    printf ("RANDOM MATRIX %iX%i\n\n", m->n, m->m);
    dump_matrix(m);
    printf ("\n\n");    
    printf(" Indeces for RANDOM MATRIX:\n");
    for (i=0;i<model->r;i++) {
      if (model->rand[i] > c) {
	c = model->rand[i];
      }
      printf(" T%i => [ %2i - %2i : ", model->rand[i], (*rindeces)[i][0], (*rindeces)[i][1]);
      for (j=0;j<(*rindeces)[i][1]-(*rindeces)[i][0];j++) {
	printf(" %2i ", model->rzeros[i][j] + 1);
      }
      printf ("]\n");
    }
    cc = model->r;
    for (i=0;i<model->i;i++) {
      c++;
      if (model->randinter[i] > 0) {
	printf(" T%i => [ %2i - %2i : ", c+i, (*rindeces)[cc][0], (*rindeces)[cc][1]);
	for (j=0;j<(*rindeces)[cc][1]-(*rindeces)[cc][0];j++) {
	  printf(" %2i ", model->rzeros[cc][j] + 1);
	}
	printf ("]\n");
	cc++;
      }
    }
    printf("\n\n");
  }

  return m;

}

MATRIX *
set_fixed_matrix(FACTOR *factor,
		 COVARIATE *covariate,
		 MODEL *model,
		 int ***findeces)
{

  int i, j, c, idx;
  int id, ty;
  MATRIX *m, *m1, *m2;

  int cc = 0;

  // Calculate the size of the findeces array
  cc = 0;
  for (i=0;i<model->f;i++) {
    if (model->infactorrand[i] == 0) {
      cc++;
    }    
  }
  for (i=0;i<model->i;i++) {
    if (model->randinter[i] == 0) {
      cc++;
    }  
  }

  *findeces = ARRAY_ALLOC(cc, int *);
  model->fzeros = ARRAY_ALLOC(cc, int *);
  model->fvalid = ARRAY_ALLOC(cc, int *);
  for (i=0;i<cc;i++) {
    (*findeces)[i] = ARRAY_ALLOC(2, int); 
  }
   
  // The first the column is the means (intercept)
  m = fill(model->s-model->k, 1, 1);
  
  // idx will have the locations of the effects
  // and interactions
  idx = 1;
  
  //Reset the counter for non random effects
  cc = 0;

  // We start with the fixed effects
  for (i=0;i<model->f;i++) {
    if (model->infactorrand[model->fixed[i]] == 0) {
      (*findeces)[cc][0] = idx;
      id = model->fixed[i];
      ty = model->tyfixed[i];
      m1 = design(factor, covariate, model, id, ty);
      m = augment_cols(m, m1, 1);
      idx = m->m;
      (*findeces)[cc][1] = idx;
      if (ty>0) {
	c = 1; 
      } else {
	c = factor->classes[id];
      }
      model->fzeros[cc] = ARRAY_ALLOC(c, int);
      model->fvalid[cc] = ARRAY_ALLOC(c, int);
      for (j=0;j<c;j++) {
        model->fzeros[cc][j] = j;
        model->fvalid[cc][j] = 1;
      }
      cc++;
    }
  }  
  
  // We continue with the interactions
  for (i=0;i<model->i;i++) {
    if (model->randinter[i] == 0) {
      m1 = fill(model->s-model->k, 1, 1);
      for (j=0;j<model->numinter[i];j++) {
	id = model->inter[i][j];
	ty = model->tyinter[i][j];
	m2 = design(factor, covariate, model, id, ty);
	m1 = cross(m1, m2, 1);
      }
      c = m1->m;
      model->fzeros[cc] = ARRAY_ALLOC(c, int);
      model->fvalid[cc] = ARRAY_ALLOC(c, int);
      m1 = clean_cross(m1, &model->fzeros[cc]);
      for (j=0;j<c;j++) {
	if (model->fzeros[cc][j] < 0) {
	  model->fvalid[cc][j] = 0;
	} else {
	  model->fvalid[cc][j] = 1;	  
	}
      }
      (*findeces)[cc][0] = idx;
      m = augment_cols(m, m1, 1);
      idx = m->m;
      (*findeces)[cc][1] = idx;
      cc++;
    }
  }
   
  if (debug) {
    printf ("FIXED MATRIX %iX%i\n\n", m->n, m->m);
    dump_matrix(m);
    printf ("\n\n");    
    printf(" Indeces for FIXED MATRIX:\n");
    for (i=0;i<cc;i++) {
      printf(" T%i => [ %2i - %2i : ", i, (*findeces)[i][0], (*findeces)[i][1]-1);
      for (j=0;j<(*findeces)[i][1]-(*findeces)[i][0];j++) {
        printf(" %2i ", model->fzeros[i][j] + 1);
      }
      printf ("]\n");          
    }
    printf("\n\n");
  }
  
  return m;

}

MATRIX *
set_design_matrix(FACTOR *factor,
		  COVARIATE *covariate,
		  MODEL *model,
		  int ***indeces)
{

  int i, j, c, idx;
  int id, ty;
  MATRIX *m, *m1, *m2;

  j = model->f + model->i;
  *indeces = ARRAY_ALLOC(j, int *);
  model->zeros = ARRAY_ALLOC(j, int *);
  model->valid = ARRAY_ALLOC(j, int *);
  for (i=0;i<j;i++) {
    (*indeces)[i] = ARRAY_ALLOC(2, int); 
  }
   
  // The first the column is the means (intercept)
  m = fill(model->s-model->k, 1, 1);

  // idx will have the locations of the effects
  // and interactions
  idx = 1;

  // We start with the fixed effects
  for (i=0;i<model->f;i++) {
    (*indeces)[i][0] = idx;
    id = model->fixed[i];
    ty = model->tyfixed[i];
    m1 = design(factor, covariate, model, id, ty);
    m = augment_cols(m, m1, 1);
    idx = m->m;
    (*indeces)[i][1] = idx;
    if (ty>0) {
      c = 1; 
    } else {
      c = factor->classes[id];
    }
    model->zeros[i] = ARRAY_ALLOC(c, int);
    model->valid[i] = ARRAY_ALLOC(c, int);
    for (j=0;j<c;j++) {
      model->zeros[i][j] = j;
      model->valid[i][j] = 1;
    }
  }  

  // We continue with the interactions
  for (i=0;i<model->i;i++) {
    m1 = fill(model->s-model->k, 1, 1);
    for (j=0;j<model->numinter[i];j++) {
      id = model->inter[i][j];
      ty = model->tyinter[i][j];
      m2 = design(factor, covariate, model, id, ty);
      m1 = cross(m1, m2, 1);
    }
    c = m1->m;
    model->zeros[model->f+i] = ARRAY_ALLOC(c, int);
    model->valid[model->f+i] = ARRAY_ALLOC(c, int);
    m1 = clean_cross(m1, &model->zeros[model->f+i]);
    for (j=0;j<c;j++) {
      if (model->zeros[model->f+i][j] < 0) {
	model->valid[model->f+i][j] = 0;
      } else {
	model->valid[model->f+i][j] = 1;	  
      }
    }
    (*indeces)[model->f+i][0] = idx;
    m = augment_cols(m, m1, 1);
    idx = m->m;
    (*indeces)[model->f+i][1] = idx;
  }
   
  if (debug) {
    printf ("DESIGN MATRIX %iX%i\n\n", m->n, m->m);
    dump_matrix(m);
    printf ("\n\n");    
    printf(" Indeces for DESIGN MATRIX:\n");
    for (i=0;i<model->f+model->i;i++) {
      printf(" T%i => [ %2i - %2i : ", i, (*indeces)[i][0], (*indeces)[i][1]-1);
      if (i < model->f) {
	for (j=0;j<(*indeces)[i][1]-(*indeces)[i][0];j++) {
	  printf(" %2i ", model->zeros[i][j] + 1);
	}
      } else {
	for (j=0;j<model->maxnumlev[i-model->f];j++) {
	  printf(" %2i ", model->zeros[i][j] + 1);
	}
      }
      printf ("]\n");          
    }
    printf("\n\n");
  }
  
  return m;

}

MATRIX *
clean_design_matrix(int *rna,
		    MATRIX *des,
		    int na)
{

  int i, j;
  MATRIX *desna;

  desna = mem_allocate_matrix(des->n - na, des->m);

  for (i=0;i<des->n-na;i++) {
    for (j=0;j<des->m;j++) {
      desna->data[i][j] = des->data[rna[i]][j];
    }
  }

  if (debug) {
    printf ("DESIGN MATRIX (no NA) %iX%i\n\n", desna->n, desna->m);
    dump_matrix(desna);
    printf ("\n\n");    
  }

  return desna;

}

NMATRIX *
set_lsmeans_matrix(MODEL *model,
		   MATRIX *des,
		   int **indeces)
{

  int i, j, k, l, m, n, nn, o, p, q, idx;
  double mean, div, tmp[1000];
  MATRIX *ls, *ls1;
  NMATRIX *lsm;

  // Memory allocation
  lsm = TYPE_ALLOC(NMATRIX);

  lsm->n = model->t;
  lsm->m = ARRAY_ALLOC(lsm->n, MATRIX *);

  // Set the matrices for each of the terms
  for (i=0;i<model->t;i++) {

    n = indeces[i][1] - indeces[i][0];
    m = des->m;
    lsm->m[i] = mem_allocate_matrix(n, m);
    ls = lsm->m[i];

    // 1. The intercept is always 1
    for (j=0;j<ls->n;j++) {
      ls->data[j][0] = 1;	
    }

    // 2. Set the values for the current term
    // whether it is a fixed effect or an
    // interaction
    if (n > 1) {
      // It is a factor
      for (j=0;j<ls->n;j++) {
	idx = 0;
	for (k=indeces[i][0];k<indeces[i][1];k++) {
	  // Set to 1 if it is the given level
	  // otherwise set to 0
	  if (j == idx) {
	    ls->data[j][k] = 1;	
	  } else {
	    ls->data[j][k] = 0;	
	  }
	  idx++;
	}
      }
    } else {
      // It is a covariate
      mean = 0.0;
      for (j=0;j<des->n;j++) {
        mean += des->data[j][n];
      }
      mean /= des->n;
      // Set to the mean of the covariate
      for (j=0;j<ls->n;j++) {
        ls->data[j][n] = mean;	
      }
    }

    // 3. Set the vslues for the other terms
    for (j=0;j<model->t;j++) {
      if (j != i) {
        // j is a fixed factor
	if (j < model->f) {
	  // i is a fixed factor
	  if (i < model->f) {
	    div = indeces[j][1] - indeces[j][0];
	    for (k=0;k<ls->n;k++) {
	      for (l=indeces[j][0];l<indeces[j][1];l++) {
		ls->data[k][l] = 1 / div;	
	      }
	    }
	  // i is an interaction
	  } else {
	    // We need to iterate over the terms again
	    // to find the dependencies, that is the
	    // terms for the interaction
	    for (k=0;k<model->t;k++) {
	      for (l=0;l<model->numdep[k];l++) {
	        if (model->dep[k][l] == i && k == j) {
      	          // These are the columns for the interaction
	          n = indeces[i][1] - indeces[i][0];
	          // These are the columns for the fixed effect
	          nn = indeces[k][1] - indeces[k][0];
	          // We iterate over the levels of the fixed effects
	          // and find how many columns are contained in the
	          // interaction
		  ls1 = lsm->m[k];
	          for (o=0;o<n;o++) {
	            for (p=0;p<nn;p++) {
		      if (ls1->data[p][indeces[i][0]+o] > 0) {
			  ls->data[o][indeces[k][0]+p] = 1;
		      } else {
			  ls->data[o][indeces[k][0]+p] = 0;
		      }
	            }
	          }
	        }
	      }
	    }
	  }
	// j is an interaction
	} else {
	  // We need to iterate over the terms again
	  // to find the dependencies, that is the
	  // terms for the interaction
	  for (k=0;k<j;k++) {
	    for (l=0;l<model->numdep[k];l++) {
	      if (model->dep[k][l] == j && k == i) {
		// These are the columns for the interaction
		n = indeces[j][1] - indeces[j][0];
		// These are the columns for the fixed effect
		nn = indeces[k][1] - indeces[k][0];
		// We iterate over the levels of the fixed effects
		// and find how many columns are contained in the
		// interaction
		for (o=0;o<nn;o++) {
		  div = 0;
		  for (p=0;p<n;p++) {
		    for (q=0;q<des->n;q++) {
		      if (des->data[q][indeces[k][0]+o] > 0 &&
			  des->data[q][indeces[j][0]+p] > 0) {
			tmp[p] = 1;
			div++;
			break;
		      } else {
			tmp[p] = 0;
		      }
		    }
		  }
		  for (p=0;p<n;p++) {
		    if (tmp[p] > 0 ) {
		      ls->data[o][indeces[j][0]+p] = 1 / div;  
		    } else {
		      ls->data[o][indeces[j][0]+p] = 0;
		    }
		  }
		}
		// We need to check is we can calculate a kosher lsmeans
		if (k < model->f) {
		  for (o=0;o<nn;o++) {
		    model->valid[k][o] = 1;
		    for (p=0;p<model->termlen[j];p++) {
		      if (model->lsm[j][p][k] == o) {
			if (model->zeros[j][p] < 0) {
			  model->valid[k][o] = 0;
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }

  }

  if (debug) {
    printf ("LEAST SQUARE MEANS MATRICES\n\n");
    for (i=0;i<lsm->n;i++) {
      printf ("TERM %i LSMEANS MATRIX %iX%i\n\n", i, lsm->m[i]->n, lsm->m[i]->m);
      dump_matrix(lsm->m[i]);
      printf ("\n");    
    }
    printf("\n\n");
  }

  return lsm;

}

MATRIX *
design(FACTOR *factor,
       COVARIATE *covariate,
       MODEL *model,
       int o,
       int t)
{

  MATRIX *d;

  int i, j, c=0;

  if (t>0) {
    d = fill(model->s, 1, 0);
    for (i=0;i<d->n;i++) {
      if (model->include[i]>0) {
	d->data[c][0] = covariate->values[o][i];
	c++;
      }
    }
  } else {
    d = fill(model->s, factor->classes[o], 0);
    for (i=0;i<d->m;i++) {
      c=0;
      for (j=0;j<d->n;j++) {
	if (model->include[j]>0) {
	  if (factor->factors[o][j]==i) {
	    d->data[c][i] = 1;
	  }
	  c++;
	}
      }
    }
  }
  // We need to adjust the number of columns
  d->n = c;

  return d;

}

int
is_nested(FACTOR *factor,
	  COVARIATE *covariate,
	  MODEL *model,
	  int f)
{

  MATRIX *f1, *f2;
  MATRIX *f1t, *f2t;
  MATRIX *f2tf1, *f1tf2;
  MATRIX *res;
  int i, j, k;
  int n=0, c=0, z=0;

  f1 = design(factor, covariate, model, f, 0);
  f1t = transpose(f1,0); 

  for (i=0;i<factor->n;i++) {
    c=0;
    if (model->infactor[i]>0 && i != f && factor->classes[i]<f1->m) {
      c++;
      f2 = design(factor, covariate, model, i, 0);
      f2t = transpose(f2,0); 
    }
    if (c>0) {
      f2tf1 = multiply(f2t,f1,0);
      f1tf2 = multiply(f1t,f2,0);
      res = multiply(f2tf1, f1tf2,1);
      z=0;
      for (j=0;j<res->n;j++) {
	for (k=0;k<res->m;k++) {
	  if (res->data[j][k]==0) {
	    z++;
	  }
	}
      }
      if (z>=(res->n*res->m/2)) {
	model->nested[f][model->numnest[f]]=i;
	model->numnest[f]++;
	model->denomnest[f]/=factor->classes[i];
	n=1;
      }
      free(f2);
      free(f2t);
      free(res);
    }
  }
  free(f1);
  free(f1t);

  return n;

}

MODEL *
set_model(FACTOR *factor,
	  COVARIATE *covariate,
	  char **fixed_effs,
	  int n_fixed_effs,
	  char **interacts,
	  int n_interacts,
	  char **random_facs,
	  int n_random_facs)
{

  int i, j, k, l, t, si, cnt, ok, max;
  int col, div, times, len, idx;
  int *tmpt, *tmpty, *tmpsrt;
  int *srtint, *lenint;
  char *tok, *ctok, str[1000];
  MODEL *m;
  
  m = TYPE_ALLOC(MODEL);

  // A. We set the fixed effects. We (1)
  // check that the name of each fixed effect
  // is either a valid factor name or a valid 
  // covariate name and (2) find all the
  // factors and covariates included in the
  // model.  We also check the names of the
  // random factors


  // 1. Allocate memory
  m->infactor = ARRAY_ALLOC(factor->n, int);
  m->infactorrand = ARRAY_ALLOC(factor->n, int);
  m->infactornest = ARRAY_ALLOC(factor->n, int);


  // 2. Initialize
  for (i=0;i<factor->n;i++) {
    // No factors are included
    m->infactor[i]=0;
    m->infactorrand[i]=0;
    m->infactornest[i]=0;
  }
  if (covariate->n > 0) {
    m->incovariate = ARRAY_ALLOC(covariate->n, int);
    for (i=0;i<covariate->n;i++) {
      // No covariates are included
      m->incovariate[i]=0;
    }
  }

  // 3. For every fixed effect name we have
  // to check all the factor and all the
  // covariate names
  m->f = n_fixed_effs;
  tmpt = ARRAY_ALLOC(m->f, int);
  tmpty = ARRAY_ALLOC(m->f, int);
  for (i=0;i<n_fixed_effs;i++) {
    ok=0;
    // check the factors first
    for (j=0;j<factor->n;j++) {
      if (!strcmp(factor->fnames[j], fixed_effs[i])) {
	// .. get the factor Id
	tmpt[i]=j;
	// .. type of fixed effect
	tmpty[i]=0;
	// .. include factor in the model
	m->infactor[j]++;
	ok++;
	break;
      }
    }
    if (ok==0) {
      // then check the covariates
      if (covariate->n > 0) {
	for (j=0;j<covariate->n;j++) {
	  if (!strcmp(covariate->cnames[j], fixed_effs[i])) {
	    // .. get the covariate id
	    tmpt[i]=j;
	    // .. type of fixed effect
	    tmpty[i]=1;
	    // include covariate in the model
	    m->incovariate[j]++;
	    ok++;
	    break;
	  }
	}
      }
      if (ok==0) {
	fprintf(stderr, "Fixed effect %s not found in factor file\n", fixed_effs[i]);
	exit(1);
      }
    }
  }

  // 4. We then sort the fixed effects
  // and set them in the model
  m->fixed = ARRAY_ALLOC(m->f, int);
  m->tyfixed = ARRAY_ALLOC(m->f, int);
  tmpsrt = ARRAY_ALLOC(m->f, int);
  sort_by_indexi(m->f, tmpt, tmpsrt);
  for (i=0;i<m->f;i++) {
    si=tmpsrt[i];
    m->fixed[i]=tmpt[si];
    m->tyfixed[i]=tmpty[si];
  }
  free( tmpt);
  free( tmpty);
  free( tmpsrt);


  // 5. We now ckeck that (1) the names of
  // the random effects are valid factor names
  // and (2) they are included in the model
  m->r = n_random_facs;
  m->rand = ARRAY_ALLOC(m->r, int);
  for (i=0;i<n_random_facs;i++) {
    ok=0;
    for (j=0;j<factor->n;j++) {
      // check the factor name
      if (!strcmp(factor->fnames[j], random_facs[i])) {
	m->rand[i] = j;
	m->infactorrand[j]++;
	ok++;
	break;
      }
    }
    if (ok == 0) {
      fprintf(stderr, "Random factor %s not found in factor file\n", random_facs[i]);
      exit(1);
    }
  }


  // B. Next, we set the interactions. We need to
  // check (1) if the terms in each interaction
  // has a valid name, either a factor or a 
  // covariate valid name, (2) include the
  // term in the model and (3) find if we have
  // enough DOF to calculate the interaction 

  // 1. We start by sorting the interactions
  // by the number of terms 
  srtint = ARRAY_ALLOC(n_interacts, int);
  lenint = ARRAY_ALLOC(n_interacts, int);
  for (i=0;i<n_interacts;i++) {
    // find number the terms in the interaction
    cnt=0;
    ctok = strdup(interacts[i]);
    tok = strtok (ctok, "*");
    while (tok != NULL) {     
      cnt++;
      tok = strtok (NULL, "*");
    }
    free(ctok);
    lenint[i]=cnt;
  }
  sort_by_indexi(n_interacts, lenint, srtint);

  // 2. Allocate memory
  m->i = n_interacts;
  m->numinter = ARRAY_ALLOC(m->i, int);
  m->maxnumlev = ARRAY_ALLOC(m->i, int);
  m->inter = ARRAY_ALLOC(m->i, int *);
  m->tyinter = ARRAY_ALLOC(m->i, int *);
  m->randinter = ARRAY_ALLOC(m->i, int);

  // 3. Check the interactions
  for (i=0;i<n_interacts;i++) {
    si=srtint[i];
    m->randinter[i] = 0;
    m->numinter[i] = lenint[si];
    max = 1;
    tmpt = ARRAY_ALLOC(1000, int);
    tmpty = ARRAY_ALLOC(1000, int);
    // check the names of all term
    cnt=0; // term number
    tok = strtok (interacts[si], "*");
    while (tok != NULL) {     
      ok=0;
      // check factor names first ..
      for (j=0;j<factor->n;j++) {
	if (!strcmp(factor->fnames[j], tok)) {
	  // .. get the factor id
	  tmpt[cnt]=j;
	  // .. type of term
	  tmpty[cnt]=0;
	  // .. include factor in the model
	  m->infactor[j]++;
	  // .. check for random factors
	  if (m->infactorrand[j]>0) {
	    m->randinter[i]++;
	  }
	  // .. get the number of levels
	  max *= factor->classes[j];
	  ok++;
	  break;
	}
      }
      if (ok == 0) {
	// now the covariate names
	if (covariate->n > 0) {
	  for (j=0;j<covariate->n;j++) {
	    if (!strcmp(covariate->cnames[j], tok)) {
	      // get the covariate id
	      tmpt[cnt]=j;
	      // type of term 
	      tmpty[cnt]=1;
	      // include the covariate in the model
	      m->incovariate[j]++;
	      ok++;
	      // .. get the number of levels which is 1
	      max *= 1;
	      break;
	    }
	  }
	}
	if (ok == 0) {
	  fprintf(stderr, "Term %s not found in factor or covariate names\n", tok);
	  exit(1);
	}
      }
      cnt++;
      tok = strtok (NULL, "*");
    }
    // Set the maximum number of levels in the interaction
    m->maxnumlev[i] = max;
    // We then sort the interaction
    // and set them in the model
    m->inter[i] = ARRAY_ALLOC(cnt, int);
    m->tyinter[i] = ARRAY_ALLOC(cnt, int);
    tmpsrt = ARRAY_ALLOC(cnt, int);
    sort_by_indexi(cnt, tmpt, tmpsrt);
    for (j=0;j<cnt;j++) {
      si=tmpsrt[j];
      m->inter[i][j]=tmpt[si];
      m->tyinter[i][j]=tmpty[si];
    }
    free(tmpt);
    free(tmpty);
    free(tmpsrt);
  }
  free(lenint);
  free(srtint);

  // C. We define the number of samples included in the
  // model based on the included factors and covariates
  // (those that are not 'Unassigned')


  // 1.  Allocate memory and set number of samples.
  // We start by inclding all samples
  m->s = factor->n > 0 ? factor->s : covariate->s;
  m->include = ARRAY_ALLOC(m->s, int);
  for (i=0;i<m->s;i++) {
    m->include[i] = 1;
  }


  // 2. We begin with the covariates. We exclude any
  // samples with missing data for a covariate.
  for (i=0;i<covariate->n;i++) {
    if (m->incovariate[i]>0) {
      for (j=0;j<m->s;j++) {
	if (isnan(covariate->values[i][j])) {
	  m->include[j] = 0;
	}
      }
    }
  }


  // 3. We continue with the factors.
  for (i=0;i<factor->n;i++) {
    if (m->infactor[i]>0) {
      for (j=0;j<m->s;j++) {
	if (factor->factors[i][j]<0) {
	  m->include[j] = 0;
	}
      }
    }
  }


  // 4. We count the number of excluded samples
  m->k = 0;
  for (i=0;i<m->s;i++) {
    if (m->include[i] < 1) {
      m->k++;
    }
  }


  // E. We check if we have nested factors so
  // we can enable back some of the required
  // interactions 
  m->n=0;
  m->infactornest = ARRAY_ALLOC(factor->n, int);
  m->numnest = ARRAY_ALLOC(factor->n, int);
  m->denomnest = ARRAY_ALLOC(factor->n, int);
  m->nested = ARRAY_ALLOC(factor->n, int *);
  for (i=0;i<factor->n;i++) {
    m->infactornest[i]=0;
    m->numnest[i]=0;
    m->denomnest[i]=factor->classes[i];
    m->nested[i] = ARRAY_ALLOC(factor->n, int);
    if (m->infactor[i]) {
      if (is_nested(factor, covariate, m, i)) {
	m->infactornest[i]++;
	m->n++;
      }
    }
  }


  // F. We finally check for dependencies
  // for the fixed effects and the 
  // interactions.


  // 1. Total number of terms in the model
  m->t = m->f + m->i;


  // 2. Allocate memory
  m->numdep = ARRAY_ALLOC(m->t, int);
  m->randdep = ARRAY_ALLOC(m->t, int);
  m->dep = ARRAY_ALLOC(m->t, int *);
  for (i=0;i<m->t;i++) {
    m->dep[i] = ARRAY_ALLOC(m->i, int);
  }


  // 3. Start with the fixed effects
  for (i=0;i<m->f;i++) {
    m->numdep[i]=0;
    for (j=0;j<m->i;j++) {
      for (k=0;k<m->numinter[j];k++) {
	if (m->fixed[i] == m->inter[j][k] &&
	    m->tyfixed[i] == m->tyinter[j][k]) {
	  m->dep[i][m->numdep[i]]=m->f+j;
	  m->numdep[i]++;
	}
      }
    }
    if (m->infactorrand[m->fixed[i]]>0) {
      m->randdep[i] = 1;
    } else {
      m->randdep[i] = 0;
    }
  }

  // 4. We continue with the interactions
  for (i=0;i<m->i;i++) {
    m->numdep[m->f+i]=0;
    for (j=i+1;j<m->i;j++) {
      if (m->numinter[j]>m->numinter[i]) {
	ok=0;
	for (k=0;k<m->numinter[i];k++) {
	  for (l=0;l<m->numinter[j];l++) {
	    if (m->inter[i][k] == m->inter[j][l] &&
		m->tyinter[i][k] == m->tyinter[j][l]) {
	      ok++;
	    }
	  }
	}
	if (ok==m->numinter[i]) {
	  m->dep[m->f+i][m->numdep[m->f+i]]=m->f+j;
	  m->numdep[m->f+i]++;
	}
      }
    }
    if (m->randinter[i]>0) {
      m->randdep[m->f+i] = 1;
    } else {
      m->randdep[m->f+i] = 0;
    }
  }


  // G. We set up the levels for the LSMEANS

  // 1. Allocate memory
  m->factors = ARRAY_ALLOC(m->t, char *);
  m->levels = ARRAY_ALLOC(m->t, char **);
  m->lsm = ARRAY_ALLOC(m->t, int **);
  m->termlen = ARRAY_ALLOC(m->t, int);

  // 2. Start with the fixed effects
  for (i=0;i<m->f;i++) {
    if (m->tyfixed[i]>0) {
      col = 1;
    } else {
      col = factor->classes[m->fixed[i]];
    }
    // a. The names of the factors
    if (col > 1) {
      // factor name
      len = strlen(factor->fnames[m->fixed[i]]);
      m->factors[i] = ARRAY_ALLOC((len+1), char);
      strncpy(m->factors[i], factor->fnames[m->fixed[i]], len);
      m->factors[i][len] = '\0';
    } else {
      // covariate name
      len = strlen(covariate->cnames[m->fixed[i]]);
      m->factors[i] = ARRAY_ALLOC((len+1), char);
      strncpy(m->factors[i], covariate->cnames[m->fixed[i]], len);
      m->factors[i][len] = '\0';
    }
    m->termlen[i] = col;
    m->levels[i] = ARRAY_ALLOC(col, char *);
    m->lsm[i] = ARRAY_ALLOC(col, int *);
    for (j=0;j<col;j++) {
      // b. The name of the levels
      if (col > 1) {
	// level names
	len = strlen(factor->cnames[m->fixed[i]][j]);
	m->levels[i][j] = ARRAY_ALLOC((len+1), char);
	strncpy(m->levels[i][j], factor->cnames[m->fixed[i]][j], len);
	m->levels[i][j][len] = '\0';
      } else {
	// covariate names
	len = strlen(covariate->cnames[m->fixed[i]]);
	m->levels[i][j] = ARRAY_ALLOC((len+1), char);
	strncpy(m->levels[i][j], covariate->cnames[m->fixed[i]], len);
	m->levels[i][j][len] = '\0';
      }
      // c. The Ids for the factors / levels
      m->lsm[i][j] = ARRAY_ALLOC(1, int);
      m->lsm[i][j][0] = j;
    }
  }

  // 3. We continue with the interactions
  for (i=0;i<m->i;i++) {
    // a. The names of the factors
    strcpy(str, "");
    for (j=0;j<m->numinter[i];j++) {
      if (m->tyinter[i][j]>0) {
	strcat(str, covariate->cnames[m->inter[i][j]]); 
      } else {
	strcat(str, factor->fnames[m->inter[i][j]]);
      }
      strcat(str, "*");
    }
    len = strlen(str)-1;
    m->factors[m->f+i] = ARRAY_ALLOC((len+1), char);
    strncpy(m->factors[m->f+i], str, len);
    m->factors[m->f+i][len] = '\0';
    // b. More memory allocation
    m->termlen[m->f+i] = m->maxnumlev[i];
    m->levels[m->f+i] = ARRAY_ALLOC(m->maxnumlev[i], char *);
    m->lsm[m->f+i] = ARRAY_ALLOC(m->maxnumlev[i], int *);
    for (j=0;j<m->maxnumlev[i];j++) {
      m->lsm[m->f+i][j] = ARRAY_ALLOC(m->numinter[i], int);	
    }
    // c. Load the Ids
    div = 1;
    for (j=0;j<m->numinter[i];j++) {
      if (m->tyinter[i][j]>0) {
	div *= 1;
        l = 1;
      } else {
        div *= factor->classes[m->inter[i][j]];
	l = factor->classes[m->inter[i][j]];
      }
      times = m->maxnumlev[i] / div;
      idx = 0;
      while (idx < m->maxnumlev[i]) {
	for (k=0;k<l;k++) {   
	  for (t=0;t<times;t++) {      
	    m->lsm[m->f+i][idx][j] = k;
	    idx++;
	  }
	}
      }
    }
    // d. Load the level combinations
    for (j=0;j<m->maxnumlev[i];j++) {
      strcpy(str, "");
      for (k=0;k<m->numinter[i];k++) {
        l = m->lsm[m->f+i][j][k];
        if (m->tyinter[i][k]>0) {
          strcat(str, covariate->cnames[m->inter[i][k]]);
        } else {
          strcat(str, factor->cnames[m->inter[i][k]][l]);
        }
        strcat(str, "*");
      }
      len = strlen(str)-1;
      m->levels[m->f+i][j] = ARRAY_ALLOC((len+1), char);
      strncpy(m->levels[m->f+i][j], str, len);
      m->levels[m->f+i][j][len] = '\0';
    }
  }

  // We are done!
  if (debug) {
    printf ("MODEL\n\n");
    printf (" Factors included in the model\n [ ");
    for (i=0;i<factor->n;i++) {
      if (m->infactor[i]>0) {
	printf("%3i ", i);
      }
    }
    printf (" ]\n\n");
    printf (" Factor names\n");
    for (i=0;i<factor->n;i++) {
      if (m->infactor[i]>0) {
	printf("%3i: %s\n", i, factor->fnames[i]);
      }
    }
    printf ("\n");
    if (covariate->n > 0) {
      printf (" Covariates included in the model\n [ ");
      for (i=0;i<covariate->n;i++) {
	if (m->incovariate[i]>0) {
	  printf("%3i ", i);
	}
      }
      printf (" ]\n\n");
      printf (" Covariate names\n");
      for (i=0;i<covariate->n;i++) {
	if (m->incovariate[i]>0) {
	  printf("%3i: %s\n", i, covariate->cnames[i]);
	}
      }
      printf ("\n");
    }
    printf (" Fixed Effects\n");
    printf (" Id(s) for fixed effects => ");
    print_array_i(m->fixed, m->f);
    printf (" Type of fixed effects   => ");
    print_array_i(m->tyfixed, m->f);
    printf ("\n");
    if (m->r > 0) {
      printf (" Random Effects\n");
      printf (" Id(s) for random effects => ");
      print_array_i(m->rand, m->r);
      printf ("\n");
    }
    if (m->i > 0) {
      printf (" Interactions\n");
      for (i=0;i<m->i;i++) {
	printf (" Ids for terms in I%-2i => ", i);
	printf (" ");
	print_array_i(m->inter[i], m->numinter[i]);	
	printf (" Type of terms in I%-2i => ", i);
	printf (" ");
	print_array_i(m->tyinter[i], m->numinter[i]);	
      }
      printf ("\n");
      printf (" Dependencies on terms\n");
      for (i=0;i<m->t;i++) {
	printf (" T%-2i => ", i); 
	if (m->numdep[i] > 0) {
	  printf (" ");
	  print_array_i(m->dep[i], m->numdep[i]);	
	} else {
	  printf (" [ None ]\n");	  
	}
      }
      printf ("\n");
    }
    printf (" Term Columns\n");
    for (i=0;i<m->t;i++) {
      printf(" T%i:\n", i); 
      for (j=0;j<m->termlen[i];j++) {
	printf("    %2i => %s\n", j+1, m->levels[i][j]); 
      }
    }
    if (m->n > 0) {
      printf (" Nested Factors\n");
      for (i=0;i<factor->n;i++) {
	if (m->infactornest[i]>0) {
	  printf (" F%-2i nested in Factor(s) ", i); 
	  print_array_i(m->nested[i], m->numnest[i]);	
	  printf (" F%-2i denominator %i\n", i, m->denomnest[i]); 
	}
      }
      printf ("\n");
    }
    printf ("\n");
    printf (" Number of Samples skipped: %i\n", m->k);
    printf ("\n");
    printf (" Samples included\n");
    printf (" ");
    print_array_i(m->include, m->s);
    printf ("\n");
    printf ("\n");
  }

  return m;

}

int *
calculate_dof(MATRIX *des,
	      MODEL *model,
	      int **indeces,
	      int ssq) 
{

  int i, *rank, *dof;
  MATRIX *gef;
  NMATRIX *hyp;
   
  // Allocate memory
  dof = ARRAY_ALLOC(model->t, int);  

  // Set the Estimable function matrix
  gef = set_gef_matrix(des, &rank);

  // Type of sum of squares
  if (ssq == 2) {
    hyp = set_ef3_matrix(model, gef, rank, indeces, 1);
  } else {
    hyp = set_ef3_matrix(model, gef, rank, indeces, 1);
  }

  for (i=0;i<model->t;i++) {
    dof[i] = hyp->m[i]->n;
  }

  free(rank);
  free(gef);
  free(hyp);

  return dof;

}

int *
get_contained_dof(MATRIX *des,
		  MODEL *model,
		  int **indeces,
		  int *ndf) 
{

  int i, j, d, r, cr;
  int rtmp, rdof, *dof;
  
  // Calculate the degrees of freedom 
  dof = ARRAY_ALLOC(model->t, int);

  // Rank for the design matrix
  r = rank_matrix(des);

  // Calculation
  for (i = 0; i < model->t; i++) {
    if (model->randdep[i] == 0) {
      if (model->numdep[i]>0) {
	rdof = r;
	for (j = 0; j < model->numdep[i]; j++) {
	  cr = 0;
	  d =  model->dep[i][j];
	  if (model->randdep[d]>0) {
	    cr++;
	    rtmp = ndf[d];
	    if (rtmp < rdof && rtmp > 0) {
	      rdof = rtmp;
	    }
	  }
	}
	if (cr>0) {
	  dof[i] = rdof;
	} else {
          dof[i] = des->n - r;
	}
      } else {
	dof[i] = des->n - r;
      }
    } else {
      dof[i] = -1;
    }
  }

  return dof;

}

double
get_sum_of_squares(MATRIX *H,
		   MATRIX *C,
		   MATRIX *b)
{

  double ss;
  MATRIX *Hb, *Hbt, *HbtC, *HbtCHb;

  Hb = multiply(H,b,0);
  Hbt = transpose(Hb,0);
  HbtC = multiply(Hbt, C, 0);
  HbtCHb = multiply(HbtC, Hb, 1);
  ss = HbtCHb->data[0][0];
  free(Hbt);
  free(HbtCHb);

  return ss;

}

MATRIX *
get_least_square_means(MATRIX *L,
		       MATRIX *B)
{

  MATRIX *lsm;

  lsm = multiply(L,B,0);

  return lsm;

}

void
parameter_estimates(VECTOR *Y,
		    MATRIX *des,
		    MATRIX *X,
		    MATRIX **b,
		    double *sse)
{

  MATRIX *y, *yt, *yty, *yhat, *ytyhat;

  // X is inv(X'X) * X'

  y = vector_to_matrix(Y,0);
  yt = transpose(y,0);
  yty = multiply(yt,y,0);
  *b = multiply(X,y,0);

  // the expected values
  yhat = multiply(des,*b,0);
  ytyhat = multiply(yt,yhat,1);
  *sse = yty->data[0][0]-ytyhat->data[0][0];
  free(y);
  free(yty);
  free(ytyhat);

}

void
satterthwaite_approx (MODEL *model,
		      double **temp,
		      double *df,
		      double *f,
		      double *ms,
		      int idx)
{

  int i, ii;
  double num=0, denom=0;

  // df are in temp[x][1]
  // ms are in temp[x][2]
  
  for (i=0;i<model->numdep[idx];i++) {
    ii = model->dep[idx][i];
    num += temp[ii][2];
    denom += pow(temp[ii][2],2)/temp[ii][1];
  }
  num = pow(num, 2);
  *df = num/denom;
  *f = temp[idx][2]/num;
  *ms = temp[idx][2];

  //ii = model->dep[idx][model->numdep[idx]-1];
  //*df = temp[ii][1];
  //*f = temp[idx][2]/temp[ii][2];
  //*ms = temp[ii][2];

}

void
glm(VECTOR *Y,
    MATRIX *design,
    MATRIX *X,
    NMATRIX *H,
    NMATRIX *C,
    NMATRIX *lsm,
    RESULT *res,
    MODEL *model,
    int idx,
    char *reference,
    int lsmo,
    int islog)
{

  int i, j, k, l, n, len, index;
  double ss, tss, ms, msm, df, sse, mse, dfe, f, pval, dft;
  double r2, ar2, aic, rv;
  double **temp;
  VECTOR *vr;
  MATRIX *b, *lsmres;
  NMATRIX *CONT,*cest,*cerr,*tstat,*fstat,*cpvalue;
  char str[1000];
  // Memory allocation
  // Estimates of contrasts
  cest = TYPE_ALLOC(NMATRIX);
  cest->n = H->n;
  cest->m = ARRAY_ALLOC(cest->n, MATRIX *);
  // St.error of contrasts
  cerr = TYPE_ALLOC(NMATRIX);
  cerr->n = H->n;
  cerr->m = ARRAY_ALLOC(cerr->n, MATRIX *);
  // T statistics of contrasts 
  tstat = TYPE_ALLOC(NMATRIX);
  tstat->n = H->n;
  tstat->m = ARRAY_ALLOC(tstat->n, MATRIX *);
  // F statistics of contrasts
  fstat = TYPE_ALLOC(NMATRIX);
  fstat->n = H->n;
  fstat->m = ARRAY_ALLOC(fstat->n, MATRIX *);
  // P -value of contrasts
  cpvalue = TYPE_ALLOC(NMATRIX);
  cpvalue->n = H->n;
  cpvalue->m = ARRAY_ALLOC(cpvalue->n, MATRIX *);

  // To store results
  vr = res->vectors[idx];

  // the parameter estimates
  n = Y->n;
  parameter_estimates(Y, design, X, &b, &sse);

  // We iterate over the hypothesis
  // to find the total DOF and the error
  // DOF
  df = 0;
  for (i=0;i<H->n-1;i++) {
    df += H->m[i]->n;
  }
  dft = n - 1;
  dfe = dft - df;

  // Calculate glm for the model
  index = H->n-1;
  df = dft - dfe;
  ss = get_sum_of_squares(H->m[H->n-1],
			  C->m[C->n-1],
			  b);

  // The error terms
  tss = sse + ss;
  mse = sse / dfe;
  ms = msm = ss/df;
  f = fabs(ms/mse);
  pval = betai(0.5*dfe,
	       0.5*df,
	       dfe/(dfe+df*f));
  vr->data[index++] = pval;
  r2 = ss/(tss);
  vr->data[index++] = r2;
  ar2 = MAX(0,1 - ((dft/(dft-df))*(1-r2)));
  vr->data[index++] = ar2;
  aic = n * log(sse/n) + (2*(b->n+1));
  vr->data[index] = aic;
  
  // Contrasts 
  CONT = set_contrast_matrices(design, H, 1);
  
  for (i=0;i<H->n-1;i++) {
    cest->m[i] = multiply(H->m[i],b,0);
    cerr->m[i] = multiply_number(CONT->m[i],mse,0);
    cerr->m[i] = extract_diag(cerr->m[i]);
    cerr->m[i] = matrix_power(cerr->m[i],0.5,0);
    tstat->m[i] = matrix_element_manipulation(cest->m[i],cerr->m[i],2);
    fstat->m[i] = matrix_element_manipulation(tstat->m[i],tstat->m[i],1);
    cpvalue->m[i] = fill(H->m[i]->n,1,0);
    for (j=0;j<H->m[i]->n;j++) {
      cpvalue->m[i]->data[j][0] = betai(0.5*dfe,
      		                        0.5,
                                        dfe/(dfe+fstat->m[i]->data[j][0]));
    }
    if (verbose || debug) {

      printf("Contrasts estimate\n");
      dump_matrix(cest->m[i]);  
      printf("Contrasts st.error\n");
      dump_matrix(cerr->m[i]); 
      printf("Contrasts t-test\n");
      dump_matrix(tstat->m[i]);
      printf("Contrasts F-test\n");
      dump_matrix(fstat->m[i]);
      printf("Contrasts p-value\n");
      dump_matrix(cpvalue->m[i]);  
    
    }

  }
  
  // Number of LSM means
  k = 0;
  for (i=0;i<lsm->n;i++) {
    for (j=0;j<model->termlen[i];j++) {
      k++;
    }
  }

  // Calculate the lsmeans and contrasts
  l = 0;
  if (reference || lsmo > 0) {
    index = model->f + model->i + 5;    
    for (i=0;i<lsm->n;i++) {
      lsmres = get_least_square_means(lsm->m[i], b);
      for (j=0;j<model->termlen[i];j++) {
	strcpy(str, model->factors[i]);
	strcat(str, ":");
	strcat(str, model->levels[i][j]);
	len = strlen(str);
	str[len]='\0';
	// Test if we have samples in the corresponding level
        if (model->zeros[i][j] > -1) {
	  // Test if we can calculate the lsmeans
	  if (model->valid[i][j] < 1) {
	    if (reference) {
	      if (!strcmp(str, reference)) {
		rv = MD;
		vr->data[index] = MD;
		l++;
	      } else {
		if (lsmo > 0 && reference) {
		  vr->data[k+index] = MD;
		  vr->data[l+index] = MD;
		} else {
		  vr->data[index] = MD;
		}
		index++;
	      }
	    } else {
	      vr->data[index++] = MD;	
	    }
	  } else {
	    if (reference) {
	      if (!strcmp(str, reference)) {
	        rv = lsmres->data[model->zeros[i][j]][0];
		vr->data[index] = lsmres->data[model->zeros[i][j]][0];
		l++;
	      } else {
		if (lsmo > 0 && reference) {
		  vr->data[k+index] = lsmres->data[model->zeros[i][j]][0];
		  vr->data[l+index] = lsmres->data[model->zeros[i][j]][0];
		} else {
		  vr->data[index] = lsmres->data[model->zeros[i][j]][0];;
		}
		index++;
	      }
	    } else {
	      vr->data[index++] = lsmres->data[model->zeros[i][j]][0];	
	    }
	  }	  
	}
      }
    }
    free(lsmres);
    if (reference) {
      j = model->f + model->i + 5;
      if (lsmo > 0) {
	j += k;
	index += k;
      }
      for (i=j;i<index;i++) {
	if (islog) {
	  vr->data[i] = vr->data[i] - rv;
	} else {
	  if (vr->data[i] > rv) {
	    vr->data[i] = vr->data[i] / rv;
	  } else {
	    vr->data[i] = (rv / vr->data[i]) * -1;
	  }
	}
      }
    }
  }

  if (verbose || debug) {
    if (reference) {
      printf("\nLEAST SQUARE MEANS (CONTRASTS)\n");
    } else {
      printf("\nLEAST SQUARE MEANS\n");
    }
    index = model->f + model->i + 5;
    for (i=0;i<lsm->n;i++) {
      lsmres = get_least_square_means(lsm->m[i], b);
      printf(" %s\n", model->factors[i]);
      for (j=0;j<model->termlen[i];j++) {
	strcpy(str, model->factors[i]);
	strcat(str, ":");
	strcat(str, model->levels[i][j]);
	len = strlen(str);
	str[len]='\0';
	if (model->zeros[i][j] > -1) {
	  // Test if we can calculate the lsmeans
	  if (model->valid[i][j] < 1) {
	    printf("%10s ", "?");	    
	    if (reference && strcmp(str, reference)) {
	      printf("(%10s) ", "?");
	    } else if (reference && !strcmp(str, reference)) {
	      printf("         ");
	    }
	  } else {
	    printf("%10.4g ", lsmres->data[model->zeros[i][j]][0]);
	    if (reference && strcmp(str, reference)) {
	      printf("(%6.2f) ", vr->data[index++]);
	    } else if (reference && !strcmp(str, reference)) {
	      printf("(%6.2f) ", 1.0);
	    }
	  }
	  printf("<= %s\n", model->levels[i][j]);
	}
      }
    }
    printf("\n");
    free(lsmres);
  }

  // Calculate glm for the terms in the model
  // and store the results in a temporary array
  index = 0;
  temp = ARRAY_ALLOC(H->n-1, double *);
  for (i=0;i<H->n-1;i++) {
    temp[i] = ARRAY_ALLOC(5, double);
    temp[i][0] = get_sum_of_squares(H->m[i],
				    C->m[i],
				    b);
    temp[i][1] = H->m[i]->n;
    temp[i][2] = temp[i][0]/temp[i][1];
    temp[i][3] = fabs(temp[i][2]/mse);
    temp[i][4] = betai(0.5*dfe,
		       0.5*temp[i][1],
		       dfe/(dfe+temp[i][1]*temp[i][3]));
    vr->data[index++] = temp[i][4];
  }

  if (verbose || debug) {
    printf("ANALYSIS OF VARIANCE                 SS     df           MS         F    Pr(>F)\n");
    printf ("%-24s %14.2f %6.0f %12.2f %9.3f  %7.4g", res->attributes[H->n-1], ss, df, ms, f, pval);      
    if (pval < 0.05) {
      printf (" *\n");
    } else {
      printf ("\n");
    }
    printf("Within (error)           %14.2f %6.0f %12.2f\n", sse, dfe, mse);
    printf("Total                    %14.2f %6.0f\n", tss, dft);
    printf("\nR2 = %6.4f     Adj.R2 = %6.4f   AIC = %6.4f\n", r2, ar2, aic);
    printf ("\nSOURCE OF VARIATION                  SS     df           MS         F    Pr(>F)\n");
    for (i=0;i<H->n-1;i++) {
      printf ("%-24s %14.2f %6.0f %12.2f %9.3f  %7.4g", res->attributes[i], temp[i][0], temp[i][1], temp[i][2], temp[i][3], temp[i][4]);      
      if (temp[i][4] < 0.05) {
	printf (" *\n");
      } else {
	printf ("\n");
      }
    }
    printf ("\n");
  }

  if (model->r > 0) {
    index=0;
    for (i=0;i<H->n-1;i++) {
      if (model->randdep[i]>0) {
	if (model->numdep[i]==0) {
	  temp[i][1] = df = dfe;
	  temp[i][3] = f = temp[i][2]/mse;
	} else {
	  satterthwaite_approx(model, temp, &df, &f, &ms, i);
	  temp[i][1] = df;
	  temp[i][2] = ms;
	  temp[i][3] = f;
	  temp[i][4] = pval = betai(0.5*dfe,
				    0.5*df,
				    dfe/(dfe+df*f));
	  vr->data[index++] = temp[i][4];
	}
      }
    }
    if (verbose || debug) {
      printf ("\nMODIFIED SOURCE OF VARIATION         ERR df     ERR MS      F    Pr(>F)\n");
      for (i=0;i<H->n-1;i++) {
	printf ("%-34s %6.0f %12.2f %9.3f  %7.4g", res->attributes[i], temp[i][1], temp[i][2], temp[i][3], temp[i][4]);	
	if (temp[i][4] < 0.05) {
	  printf (" *\n");
	} else {
	  printf ("\n");
	}
      }
      printf ("\n");
    }
  }

  free(temp);
  free(b);

}

void
mixed(VECTOR *y,
      MATRIX *z,
      int **rindeces,
      MATRIX *f,
      MATRIX *R,
      MATRIX *G,
      MATRIX **hess,
      MATRIX *grad,
      NMATRIX *G_deriv,
      NMATRIX *R_deriv,
      MATRIX *design,
      int **indeces,
      MATRIX *X,
      NMATRIX *H,
      int *ndf,
      NMATRIX *C,
      NMATRIX *lsm,
      RESULT *res,
      MODEL *model,
      enum METHOD method,
      int idx,
      char *reference,
      int lsmo,
      int islog)
{
 
  int i, *rank, rr, niter, *ddf, *temp2;
  double det, sse, reml, ml, mvq, minm;
  MATRIX *parameters;
  MATRIX *Y, *b;
  MATRIX *iR;
  MATRIX *iW, *W0, *Wa, *W, *bu, *varb;

  int j, k, n, r, index;
  double df, dfs, mse, dfe, pval, dft;
  double **temp, *temp1;
  VECTOR *vr;
  MATRIX *lsmres;
  char str[1000];

  int s, e, l, len;
  double rv;

  // To store results
  vr = res->vectors[idx];

  // Convert the vector to a matrix to simplify calculations
  Y = vector_to_matrix(y,0);
  
  // Get the parameter estimates
  n = y->n;

  // Set parameter acording to GLM if we have REML or ML
  parameter_estimates(y, f, X, &b, &sse);

  // We iterate over the hypothesis to find the total DOF and the
  // error DOF
  df = 0;
  for (i=0;i<H->n-1;i++) {
    df += H->m[i]->n;
  }
  dft = n - 1;
  dfe = dft - df;
  mse = sse / dfe;

  // Mixed Model Equations according to:
  // SAS Mixed Models. Page 743

  if (method == MIVQUE) {
    // Set the initial W Matrix
    iW = set_initial_W_matrix(f,z,Y,&bu,&varb);
 
    // Compute the initial derivatives
    niter = 1;
    compute_derivatives(f,z,R,G,Y,&*hess,&grad,iW,bu,varb,G_deriv,R_deriv,method,niter);
 
    // Get the initial parameters (MIVQUE 0)
    parameters = get_Newton_direction_matrix(&*hess, grad,0);
 
    for (i=0; i< parameters->n; i++) {  
      if (parameters->data[i][0]<0) {
	parameters->data[i][0] = 0;
      }
    }
  
  } else {

      // Set the initial W Matrix
      iW = set_initial_W_matrix(f,z,Y,&bu,&varb);

      // Compute the initial derivatives
      niter = 1;
      compute_derivatives(f,z,R,G,Y,&*hess,&grad,iW,bu,varb,G_deriv,R_deriv,method,niter);

      // Get the initial parameters (MIVQUE 0)
      parameters = get_Newton_direction_matrix(&*hess, grad,0);
      minm = min_matrix(parameters);
  
      for (i=0; i< parameters->n; i++) {  
	  if (parameters->data[i][0]<0) {
	      parameters->data[i][0] = 0;
	  }
      }
    
 
  }

  // Update G and R matrices
  update_G_matrix(&G,parameters,rindeces);
  update_R_matrix(&R,parameters);

  // We invert the r Matrix only
  // once to improve performance
  rank = ARRAY_ALLOC(R->m, int);
  iR = g2invert(R,
		R->m,
		&rank,
		&rr,
		&det,
		0);
  // Compute Likelihood Calculations   
  W0 = compute_W0_transform(f,z,iR,Y);
  Wa = compute_W0_augment(W0,iR,f,z,G);
  W = compute_W_transform(Wa,G,f,R,&bu,&varb,method,&ml,&reml,0);
  mvq = reml;
      
  if (method != MIVQUE) {
      Expected_Maxim(f,z,iR,G,R,rindeces,Y,parameters,&*hess,&grad,&bu,&varb,G_deriv,R_deriv,method,&ml,&reml);
      //Newton_Raphson(f,z,iR,G,R,rindeces,Y,parameters,&*hess,&grad,&bu,&varb,G_deriv,R_deriv,method,&ml,&reml);
  }
        
  // Calculate F for the model
  b = copy(bu,0,0,f->m,bu->m);
  index = model->t - model->r;

  // Calculate the lsmeans and contrasts
  l = 0;
  if (reference || lsmo > 0) {
    s = model->t + 3;
    e = s + model->f - model->r;
    for (i=0;i<model->i;i++) {
      if (model->randinter[i] == 0) {
        e++;
      }
    }
    index = e;
    for (i=0;i<lsm->n;i++) {
      lsmres = get_least_square_means(lsm->m[i], bu);
      for (j=0;j<model->termlen[i];j++) {
        strcpy(str, model->factors[i]);
        strcat(str, ":");
        strcat(str, model->levels[i][j]);
        len = strlen(str);
        str[len]='\0';
        // Test if we have samples in the corresponding level
        if (model->zeros[i][j] > -1) {
          // Test if we can calculate the lsmeans
          if (model->valid[i][j] < 1) {
            if (reference) {
              if (!strcmp(str, reference)) {
        	rv = MD;
        	vr->data[index] = MD;
        	l++;
              } else {
        	if (lsmo > 0 && reference) {
        	  vr->data[k+index] = MD;
        	  vr->data[l+index] = MD;
        	} else {
        	  vr->data[index] = MD;
        	}
        	index++;
              }
            } else {
              vr->data[index++] = MD;	
            }
          } else {
            if (reference) {
              if (!strcmp(str, reference)) {
                rv = lsmres->data[model->zeros[i][j]][0];
        	vr->data[index] = lsmres->data[model->zeros[i][j]][0];
        	l++;
              } else {
        	if (lsmo > 0 && reference) {
        	  vr->data[k+index] = lsmres->data[model->zeros[i][j]][0];
        	  vr->data[l+index] = lsmres->data[model->zeros[i][j]][0];
        	} else {
        	  vr->data[index] = lsmres->data[model->zeros[i][j]][0];;
        	}
        	index++;
              }
            } else {
              vr->data[index++] = lsmres->data[model->zeros[i][j]][0];	
            }
          }	  
        }
      }
    }
    free(lsmres);
    if (reference) {
      j = model->f + model->i + 5;
      if (lsmo > 0) {
        j += k;
        index += k;
      }
      for (i=j;i<index;i++) {
        if (islog) {
          vr->data[i] = vr->data[i] - rv;
        } else {
          if (vr->data[i] > rv) {
            vr->data[i] = vr->data[i] / rv;
          } else {
            vr->data[i] = (rv / vr->data[i]) * -1;
          }
        }
      }
    }
  }
  
  if (verbose || debug) {
    if (reference) {
      printf("\nLEAST SQUARE MEANS (CONTRASTS)\n");
    } else {
      printf("\nLEAST SQUARE MEANS\n");
    }
    index = model->f + model->i + 5;
    for (i=0;i<lsm->n;i++) {
      lsmres = get_least_square_means(lsm->m[i], bu);
      printf(" %s\n", model->factors[i]);
      for (j=0;j<model->termlen[i];j++) {
        strcpy(str, model->factors[i]);
        strcat(str, ":");
        strcat(str, model->levels[i][j]);
        len = strlen(str);
        str[len]='\0';
        if (model->zeros[i][j] > -1) {
          // Test if we can calculate the lsmeans
          if (model->valid[i][j] < 1) {
            printf("%10s ", "?");	    
            if (reference && strcmp(str, reference)) {
              printf("(%10s) ", "?");
            } else if (reference && !strcmp(str, reference)) {
              printf("         ");
            }
          } else {
            printf("%10.4g ", lsmres->data[model->zeros[i][j]][0]);
            if (reference && strcmp(str, reference)) {
              printf("(%6.2f) ", vr->data[index++]);
            } else if (reference && !strcmp(str, reference)) {
              printf("(%6.2f) ", 1.0);
            }
          }
          printf("<= %s\n", model->levels[i][j]);
        }
      }
    }
    printf("\n");
    free(lsmres);
  }

  index = 0;
  temp = ARRAY_ALLOC(H->n-1, double *);
  temp1 = ARRAY_ALLOC(H->n-1, double);
  temp2 = ARRAY_ALLOC(H->n-1, int);

  ddf = get_contained_dof(design,model,indeces,ndf);

  k = 0;
  for (i=0;i<model->t;i++) {
    if (model->randdep[i] == 0) {
      temp1[k] = ddf[i];
      temp2[k] = i;
      k++;
    }
  }

  // F statistics computation for fixed factors and interactions

  for (i=0;i<H->n-1;i++) {
      dfs = compute_f_statistic(H->m[i], varb, b);
      r = rank_matrix(H->m[i]);
      df = H->m[i]->n;
      dfe = temp1[i];
      vr->data[index] = betai(0.5*dfe,0.5*df, dfe/(dfe+df*dfs/r)); 
      index++;
  }

  // Variance component calculations for the random factors
  for (i=0; i< parameters->n; i++) {  
      vr->data[index++] = parameters->data[i][0];
  }
  
  // AIC calculation
  k = 0;
  for (i=0; i< parameters->n; i++) {  
    if (parameters->data[i][0]> TOL) {
      k++;
    }
  }      

  // Likelihood and AIC 
  if (method == REML) {
      vr->data[index++] = reml+2*(k);
      vr->data[index++] = reml;
  } else if (method == ML) {
      vr->data[index++] = ml+2*(k); 
      vr->data[index++] = ml;
  } else if (method == MIVQUE) {
      vr->data[index++] = mvq+2*(k); 
      vr->data[index++] = mvq;  
  }
  
  if (verbose || debug) {
      printf ("\nEffects         Num df     Den df	  F    	Pr(>F)\n");
      for (i=0;i<H->n-1;i++) {
	  dfs = compute_f_statistic(H->m[i], varb, b);
	  r = rank_matrix(H->m[i]);
	  df = H->m[i]->n;
	  dfe = temp1[i];
	  pval = betai(0.5*dfe,0.5*df, dfe/(dfe+df*dfs/r)); 
	  printf ("%-18s %3.0f %10.0f %9.3f  %9.4g", res->attributes[i], df, dfe, dfs/r, pval);	
	  if (pval < 0.05) {
	      printf (" *\n");
	  } else {
	      printf ("\n");
	  }
      }
      printf ("\n");
  }

  free(temp);
  free(b);

}

float
compute_f_statistic(MATRIX *H,
		    MATRIX *varb,
		    MATRIX *b)
{

    int *rank, rr;
    double lv, dfs;
    MATRIX *tmp1, *tmp2, *tmp3, *tmp4, *tmp5, *tmp6, *tmp7, *tmp8;

    tmp1 = multiply(H, varb, 0);
    tmp2 = transpose(H, 0);
    tmp3 = multiply(tmp1, tmp2, 0);
    rank = ARRAY_ALLOC(tmp3->m, int);  
    tmp3 = g2invert(tmp3,
		    tmp3->m,
		    &rank,
		    &rr,
		    &lv,
		    0);
    tmp4 = multiply(tmp2, tmp3, 0);
    tmp5 = multiply(tmp4, H, 0);
    tmp6 = transpose(b,0);
    tmp7 = multiply(tmp6, tmp5, 0);
    tmp8 = multiply(tmp7, b, 0);

    dfs = tmp8->data[0][0];

    free(tmp1);
    free(tmp2);
    free(tmp3);
    free(tmp4);
    free(tmp5);
    free(tmp6);
    free(tmp7);
    free(tmp8);
    free(rank);

    return dfs;

}

int 
main(int argc,
     char **argv)

{

  char *data_file, *factor_file, *covariate_file, *output_file, *reference, **variable_names;
  int  samples, variables, n_fixed_effs, n_interacts, n_random_facs, ssq, permuts, lines, variable_number, lsmo, islog;
  char **fixed_effs, **interacts, **random_facs; 
  double Q;
  enum METHOD method;
  XDATA *xdata;
  FACTOR *factor;
  COVARIATE *covariate;
  MODEL *model;
  RESULT *res;
  
  verbose = 0;

  debug = 0;

  parse_arguments(argc, argv,
		  &data_file,
		  &factor_file, &covariate_file,
		  &output_file,
		  &samples, &variables,
		  &variable_names, &variable_number,
		  &fixed_effs, &n_fixed_effs,
		  &interacts, &n_interacts,
		  &random_facs, &n_random_facs, &method,
		  &ssq, &Q, &reference, &lsmo, &islog, &permuts);
    
  if (factor_file != NULL) {

    if (samples == 0) {
      samples = count_fields(factor_file);
    }

    lines = count_lines(factor_file);
    lines--;
  
  } else {

    lines = 0;

  }

  factor = read_factor(factor_file, lines, samples);

  if (covariate_file != NULL) {

    if (samples == 0) {
      samples = count_fields(covariate_file);
    }

    lines = count_lines(covariate_file);
    lines--;
  
  } else {

    lines = 0;

  }

  covariate = read_covariate(covariate_file, lines, samples);

  model = set_model(factor, covariate,
		    fixed_effs, n_fixed_effs,
		    interacts, n_interacts,
		    random_facs, n_random_facs);

  if (reference != NULL) {

    validate_reference(model, reference);

  }

  if (variables == 0) {
    variables = count_lines(data_file);
    variables--;
  }
  
  xdata = read_xdata(data_file, samples, variables);

  res = calculate_model(xdata, factor, covariate,
			variable_names, variable_number,
			model, method, ssq, Q, reference,
			lsmo, islog, permuts);

  print_results(output_file, res, 1, model->t);
  
  free(data_file);
  free(factor_file);
  free(covariate_file);
  free(output_file);
  free(variable_names);
  free(fixed_effs);
  free(interacts);
  free(random_facs); 
  free(xdata);
  free(factor);
  free(covariate);
  free(model);
  free(res);

  exit(0);

}
