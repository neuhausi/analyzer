/*********************************************************************
 Copyright 2007 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

 $Id: nipals.c,v 1.9 2011/03/17 22:11:18 neuhausi Exp $
**********************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#define EXTERN
#include "nipals.h"

void 
parse_arguments(int argc,
                char **argv,
                char **data_file,
                char **factor_file,
                char **output_file,
                int *samples,
                int *variables,
                int *components,
                int *tolerance,
                int *all,
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

  *tolerance = 12;

  *all = 0;

  *transpose = 0;

  while ((c = getopt(argc, argv, "d:f:o:s:g:c:p:atvD")) != EOF)
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
    case 'p':
      *tolerance = atoi(optarg);
      break;
    case 'a':
      *all = 1;
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
    printf("\nUsage: nipals <options> (options in square brackets are optional)\n\n");
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
    printf(" -p [<integer>] negative power to specify tolerance (defaults to 12).\n");
    printf(" -a [<switch>]  output loadings, sum of squares, residuals and leverages.\n");
    printf(" -t [<switch>]  transpose data matrix.\n\n");
    printf("Examples\n");
    printf("\n\tnipals -d ./t/nipals.dat\n\n");
    exit(0);
  }
   
}

RESULT **
run_nipals(MATRIX *matrix,
           int components,
           int tolerance,
           int all)
{

  int i, j, k, n, m, md, max, iter;
  double TxT, Pn, c;
  double ssqX, *ssq, *ssqM, *ssqMr, *ssqMc, *T_old;
  double **T, **P, **Rs, **Rv;
  double TOLERANCE;

  RESULT **res;
  VECTOR *vr;

  // tolerance
  TOLERANCE = pow(10,tolerance * -1);

  // Size of the matrix
  n = matrix->n;
  m = matrix->m;

  // Maximum number of components
  max = components > m ? m : components;

  // Memory Allocation and matrix initialization
  T = ARRAY_ALLOC(m, double *);
  Rs = ARRAY_ALLOC(m, double *);
  for (i=0;i<m;i++) {
    T[i] = ARRAY_ALLOC(max, double);
    Rs[i] = ARRAY_ALLOC(max, double);
    for (j=0;j<max;j++) {
      T[i][j]  = 0;
      Rs[i][j] = 0;
    }
  }
  P = ARRAY_ALLOC(n, double *);
  Rv = ARRAY_ALLOC(n, double *);
  for (i=0;i<n;i++) {
    P[i] = ARRAY_ALLOC(max, double);
    Rv[i] = ARRAY_ALLOC(max, double);
    for (j=0;j<max;j++) {
      P[i][j]  = 0;
      Rv[i][j] = 0;
    }
  }
  ssq = ARRAY_ALLOC(max, double);
  for (i=0;i<max;i++) {
    ssq[i] = 0;
  }
  ssqM = ARRAY_ALLOC(n, double);
  ssqMr = ARRAY_ALLOC(n, double);
  for (i=0;i<n;i++) {
    ssqM[i] = 0;
    ssqMr[i] = 0;
  }
  ssqMc = ARRAY_ALLOC(m, double);
  T_old = ARRAY_ALLOC(m, double);
  for (i=0;i<m;i++) {
    ssqMc[i] = 0;
    T_old[i] = 0;
  }

  // Set the Result object
  res = set_res_object(matrix, components, all);

  // Scale data
  zscore(matrix->data, matrix->n, matrix->m, 1);

  // Check if we have missing data
  md = 0;
  for (i=0;i<n;i++) {
    if (matrix->md[i] > 0) {
      md += matrix->md[i];      
    }
  }

  // Calculate Sum of Squares 
  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
      if (! isnan(matrix->data[i][j])) {
         ssqM[i] += pow(matrix->data[i][j],2);
      }
    }
  }
  ssqX = 0.0;
  for (i=0;i<n;i++) {
    ssqX += ssqM[i];  
  }

  // Iterate over the components required
  for (i=0;i<max;i++) {
    iter = 0;
    for (j=0;j<m;j++) {
      if (! isnan(matrix->data[0][j])) {
        T[j][i] = matrix->data[0][j];
        T_old[j] = T[j][i] * 100;
      }
    }
    if (md > 0) {
      while ((is_continue(T_old,T,m,i) > TOLERANCE) && (iter < MAX_ITER)) {
        iter++;
        for (j=0;j<m;j++) {
          T_old[j] = T[j][i];
        }
        for (j=0;j<n;j++) {
          P[j][i] = 0;
        }
        for (j=0;j<n;j++) {
          for (k=0;k<m;k++) {
	    if (! isnan(matrix->data[j][k]) && ! isnan(T[k][i])) {            
              P[j][i] += (matrix->data[j][k] * T[k][i]);
            }
          }
        }   
        for (j=0;j<n;j++) {
          c = 0;
          for (k=0;k<m;k++) {
            if (! isnan(T[k][i])) {            
              c +=  pow(T[k][i],2);
            }
          }
          if (c) {
            P[j][i] /= c;
          }
        }
        Pn = 0;
        for (j=0;j<n;j++) {
          Pn += P[j][i] * P[j][i];
        }
        Pn = pow(Pn,0.5);
        for (j=0;j<n;j++) {
          P[j][i] /= Pn;
        }
        for (j=0;j<m;j++) {
          T[j][i] = 0;
        }
        for (j=0;j<n;j++) {
          for (k=0;k<m;k++) {
            if (! isnan(matrix->data[j][k]) && ! isnan(P[j][i])) {            
              T[k][i] += matrix->data[j][k] * P[j][i] * pow(m-1,0.5)/pow(m,0.5);
            }
          }
        }   
        for (j=0;j<m;j++) {
          c = 0;
          for (k=0;k<n;k++) {
            if (! isnan(P[k][i])) {            
              c +=  pow(P[k][i],2);
            }
          }
          if (c) {
            T[j][i] /= c;
          }
        }
      }
    } else {
      while ((is_continue(T_old,T,m,i) > TOLERANCE) && (iter < MAX_ITER)) {
        iter++;
        for (j=0;j<m;j++) {
          T_old[j] = T[j][i];
        }
        TxT = 0;
        for (j=0;j<m;j++) {
          TxT += T[j][i] * T[j][i];
        }
        for (j=0;j<n;j++) {
          P[j][i] = 0;
        }
        for (j=0;j<n;j++) {
          for (k=0;k<m;k++) {
            P[j][i] += (matrix->data[j][k] * T[k][i]) / TxT;
          }
        }   
        Pn = 0;
        for (j=0;j<n;j++) {
          Pn += P[j][i] * P[j][i];
        }
        Pn = pow(Pn,0.5);
        for (j=0;j<n;j++) {
          P[j][i] /= Pn;
        }
        for (j=0;j<m;j++) {
          T[j][i] = 0;
        }
        for (j=0;j<n;j++) {
          for (k=0;k<m;k++) {
            T[k][i] += (matrix->data[j][k] * P[j][i]) * pow(m-1,0.5)/pow(m,0.5);
          }
        }   
      }
    }
    if (iter == MAX_ITER) {
      exit_on_error("Maximum number iterations reached before convergence");    
    }
    for (j=0;j<n;j++) {
      for (k=0;k<m;k++) {
	if (! isnan(matrix->data[j][k]) && ! isnan(T[k][i]) && ! isnan( P[j][i])) {
          matrix->data[j][k] -= (T[k][i] * P[j][i]);
        }
      }
    }   
    if (all == 1) {
      // Sum of Squares for the modified matrix
      for (j=0;j<m;j++) {
        ssqMc[j] = 0;
        for (k=0;k<n;k++) {
          if (! isnan(matrix->data[k][j])) {
            ssqMc[j] += pow(matrix->data[k][j],2);
          }
        }
      }
      for (j=0;j<n;j++) {
        ssqMr[j] = 0;
        for (k=0;k<m;k++) {
          if (! isnan(matrix->data[j][k])) {
            ssqMr[j] += pow(matrix->data[j][k],2);
          }
        }
      }
      // Sum of Squares for the component
      for (j=0;j<m;j++) {
        ssq[i] += ssqMc[j]; 
      }
      ssq[i] = (ssqX - ssq[i]) / ssqX;
      // Sample Residuals
      for (j=0;j<m;j++) {
        Rs[j][i] = sqrt(ssqMc[j]);
      }
      // Variable Residuals
      for (j=0;j<n;j++) {
        Rv[j][i] = sqrt(ssqMr[j]);
      }
    }
  }

  // Set the result object

  // Scores
  for (i=0;i<m;i++) {
    vr = res[0]->vectors[i];      
    for (j=0;j<max;j++) {
      vr->data[j] = T[i][j];
    }
  }
  if (all == 1) {
    // Loadings
    for (i=0;i<n;i++) {
      vr = res[1]->vectors[i];      
      for (j=0;j<max;j++) {
        vr->data[j] = P[i][j];
      }
    }
    // Sum of Squares
    for (i=0;i<max;i++) {
      vr = res[2]->vectors[i];      
      for (j=0;j<1;j++) {
        vr->data[j] = ssq[i];
      }
    }
    // Sample Residuals
    for (i=0;i<m;i++) {
      vr = res[3]->vectors[i];      
      for (j=0;j<max;j++) {
        vr->data[j] = Rs[i][j];
      }
    }
    // Variable Residuals
    for (i=0;i<n;i++) {
      vr = res[4]->vectors[i];      
      for (j=0;j<max;j++) {
        vr->data[j] = Rv[i][j];
      }
    }
  }

  free(ssq);
  free(ssqM);
  free(ssqMr);  
  free(ssqMc);
  free(T_old);
  free(T);
  free(P);
  free(Rs);
  free(Rv);


  return res;
 
}

double
is_continue (double *T_old,
             double **T,
             int m,
             int idx)
{

  int i;
  double ssq = 0;

  for (i=0;i<m;i++) {
    if (! isnan(T_old[i])) {
      ssq += pow(T_old[i] - T[i][idx], 2);
    }
  }

  return ssq;

}

RESULT **
set_res_object(MATRIX *matrix,
               int c,
               int a)
{

  int i, j, n, s, max, nres, len;
  RESULT **res;
  VECTOR *v;
  char strc[1000];
  char str[1000];

  // The size of the matrix
  s = matrix->m;
  n = matrix->n;

  // Maximum number of components
  max = c > s ? s : c;

  // Number of results
  nres = a == 1 ? 5 : 1; 

  // Allocate results
  res = ARRAY_ALLOC(nres, RESULT *);
  for (i=0; i<nres; i++) {
    res[i] = TYPE_ALLOC(RESULT);
  }

  // Results

  for (i=0; i<nres; i++) {
    if (i == 0 || i == 3) {
      // (Samples x Components)
      // 0: Scores
      // 3: Sample Residuals
      res[i]->n = s;
      res[i]->a = max;
    } else if (i == 1 || i == 4) {
      // (Variables x Components)
      // 1: Loadings
      // 4: Variable Residuals
      res[i]->n = n;
      res[i]->a = max;
    } else if (i == 2) {          
      // (Components x 1)
      // 2: Cumulative Sum of Squares
      res[i]->n = max;
      res[i]->a = 1;
    }
    // Allocate memory
    res[i]->attributes = ARRAY_ALLOC(res[i]->a, char *);
    res[i]->variables = ARRAY_ALLOC(res[i]->n, char *);
    res[i]->vectors = ARRAY_ALLOC(res[i]->n, VECTOR *);
    for (j=0;j<res[i]->n;j++) {
      v = res[i]->vectors[j] = TYPE_ALLOC(VECTOR);
      v->data = ARRAY_ALLOC(res[i]->a, double);
      v->n = res[i]->a;
    }
    // Load the names for attributes
    if (i != 2) {
      for (j=0;j<res[i]->a;j++) {
        strcpy(str, "PC");
        sprintf(strc, "%i", j+1);
        strcat(str, strc);
        len = strlen(str);
        res[i]->attributes[j] = ARRAY_ALLOC((len+1), char);
        strncpy(res[i]->attributes[j], str, len);
        res[i]->attributes[j][len] = '\0';
      }
    } else if (i == 2) {
      strcpy(str, "Cumulative Sum of Squares");
      len = strlen(str);
      res[i]->attributes[0] = ARRAY_ALLOC((len+1), char);
      strncpy(res[i]->attributes[0], str, len);
      res[i]->attributes[0][len] = '\0';
    }
    // Load the names for variables
    if (i == 0 || i == 3) {
      for (j=0;j<res[i]->n;j++) {
        strcpy(str, matrix->cols[j]);
        len = strlen(str);
        res[i]->variables[j] = ARRAY_ALLOC((len+1), char);
        strncpy(res[i]->variables[j], str, len);
        res[i]->variables[j][len] = '\0';
      }
    } else if (i == 1 || i == 4) {
      for (j=0;j<res[i]->n;j++) {
        strcpy(str, matrix->rows[j]);
        len = strlen(str);
        res[i]->variables[j] = ARRAY_ALLOC((len+1), char);
        strncpy(res[i]->variables[j], str, len);
        res[i]->variables[j][len] = '\0';
      }
    } else if (i == 2) {          
      for (j=0;j<res[i]->n;j++) {
        strcpy(str, "PC");
        sprintf(strc, "%i", j+1);
        strcat(str, strc);
        len = strlen(str);
        res[i]->variables[j] = ARRAY_ALLOC((len+1), char);
        strncpy(res[i]->variables[j], str, len);
        res[i]->variables[j][len] = '\0';
      }
    }
  }
  
  return res;

}

char **
set_out_files(char *output_file,
              int nres)
{

  int i;
  char **outfiles;

  outfiles = ARRAY_ALLOC(nres, char *);

  if (output_file != NULL) {
    for (i=0;i<nres;i++) {
      outfiles[i] = ARRAY_ALLOC((strlen(output_file) + strlen(".xxx") + 1), char);
      strcpy(outfiles[i], output_file);
      if (i == 0) {
        strcat(outfiles[i], ".scr");
      } else if (i == 1) {
        strcat(outfiles[i], ".ldg");
      } else if (i == 2) {
        strcat(outfiles[i], ".ssq");
      } else if (i == 3) {
        strcat(outfiles[i], ".srs");
      } else if (i == 4) {
        strcat(outfiles[i], ".vrs");
      }
    }
  } else {
    for (i=0;i<nres;i++) {
      outfiles[i] = NULL;
    }
  }

  return outfiles;

}

int
main(int argc,
     char **argv)

{

    char *data_file, *factor_file, *output_file, **outfiles;
    int  samples, variables, components, tolerance, all, transpose, binary;
    int nres;
    MATRIX *matrix;
    FACTOR *factor;
    RESULT **res;

    verbose = 0;

    debug = 0;

    binary = 0;

    parse_arguments(argc, argv, &data_file, &factor_file, &output_file, 
                    &samples, &variables, &components, &tolerance,
                    &all, &transpose);

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

    res = run_nipals(matrix, components, tolerance, all);

    nres = all == 1 ? 5 : 1; 

    outfiles = set_out_files(output_file, nres);

    print_mult_results(nres, outfiles, res, 0, 0); 
    
    free(data_file);
    free(factor_file);
    free(output_file);
    free(matrix);
    free(res);

    exit(0);

}
