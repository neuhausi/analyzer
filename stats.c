/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: stats.c,v 1.26 2008/09/11 17:39:47 neuhausi Exp $
**********************************************************************/

#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <unistd.h>
#include <math.h>
#define EXTERN
#include "stats.h"

void
simple_sort(int n, double *arr)
{
  int i, ir=n-1, j, k, l=0, istack[n], jstack=0;
  double a, temp;

  for (;;) {
    if (ir-l < 8) {
      for (j = l+1; j <= ir; j++) {
	a = arr[j];
	for (i = j-1; i >= l; i--) {
	  if (arr[i] <= a) break;
	  arr[i+1] = arr[i];
	}
	arr[i+1] = a;
      }
      if (jstack == 0) break;
      ir = istack[jstack--];
      l = istack[jstack--];
    } else {
      k = (l+ir) >> 1;
      SWAP(arr[k], arr[l+1])
	if (arr[l] >  arr[ir]) {
	  SWAP(arr[l], arr[ir])
	  }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1], arr[ir])
	}
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l], arr[l+1])
	}
      i = l+1;
      j = ir;
      a = arr[l+1];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i], arr[j]);
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      jstack += 2;
      if (jstack >= n) exit_on_error("stack too small in sort.");
      if (ir-i+1 >= j-l) {
	istack[jstack] = ir;
	istack[jstack-1] = i;
	ir = j-1;
      } else {
	istack[jstack] = j-1;
	istack[jstack-1] = l;
	l = i;
      }
    }
  }

}

void
sort (int n, double *uar, double *sar)
{
  int i, ir=n-1, j, k, l=0, istack[n], jstack=0;
  double a, b, temp;
  
  for (;;) {
    if (ir-l < 8) {
      for (j = l+1; j <= ir; j++) {
	a = uar[j];
	b = sar[j];
	for (i = j-1; i >= l; i--) {
	  if (uar[i] <= a) break;
	  uar[i+1] = uar[i];
	  sar[i+1] = sar[i];
	}
	uar[i+1] = a;
	sar[i+1] = b;
      }
      if (!jstack) {
	return;
      }
      ir = istack[jstack];
      l  = istack[jstack-1];
      jstack -= 2;
    } else {
      k = (l+ir) >> 1;
      SWAP(uar[k], uar[l+1])
	SWAP(sar[k], sar[l+1])
	if (uar[l] > uar[ir]) {
	  SWAP(uar[l],uar[ir])
	    SWAP(sar[l],sar[ir])
	    }
      if (uar[l+1] > uar[ir]) {
	SWAP(uar[l+1], uar[ir])
	  SWAP(sar[l+1], sar[ir])
	  }
      if (uar[l] > uar[l+1]) {
	SWAP(uar[l], uar[l+1])
	  SWAP(sar[l], sar[l+1])
	  }
      i = l+1;
      j = ir;
      a = uar[l+1];
      b = sar[l+1];
      for (;;) {
	do i++; while (uar[i] < a);
	do j--; while (uar[j] > a);
	if (j < i) break;
	SWAP(uar[i], uar[j])
	  SWAP(sar[i], sar[j])
	  }
      uar[l+1] = uar[j];
      uar[j] = a;
      sar[l+1] = sar[j];
      sar[j] = b;
      jstack += 2;
      if (jstack >= n) exit_on_error("stack too small in sort");
      if (ir-i+1 >= j-l) {
	istack[jstack] = ir;
	istack[jstack-1] = i;
	ir = j-1;
      } else {
	istack[jstack] = j-1;
	istack[jstack-1] = l;
	l = i;
      }
    }
  }

}

void
sorti (int n, int *uar, int *sar)
{
  int i, ir=n-1, j, k, l=0, istack[n], jstack=0;
  int a, b, temp;
  
  for (;;) {
    if (ir-l < 8) {
      for (j = l+1; j <= ir; j++) {
	a = uar[j];
	b = sar[j];
	for (i = j-1; i >= l; i--) {
	  if (uar[i] <= a) break;
	  uar[i+1] = uar[i];
	  sar[i+1] = sar[i];
	}
	uar[i+1] = a;
	sar[i+1] = b;	
      }
      if (!jstack) {
	return;
      }
      ir = istack[jstack];
      l  = istack[jstack-1];
      jstack -= 2;
    } else {
      k = (l+ir) >> 1;
      SWAP(uar[k], uar[l+1])
	SWAP(sar[k], sar[l+1])
	if (uar[l] > uar[ir]) {
	  SWAP(uar[l],uar[ir])
	    SWAP(sar[l],sar[ir])
	    }
      if (uar[l+1] > uar[ir]) {
	SWAP(uar[l+1], uar[ir])
	  SWAP(sar[l+1], sar[ir])
	  }
      if (uar[l] > uar[l+1]) {
	SWAP(uar[l], uar[l+1])
	  SWAP(sar[l], sar[l+1])
	  }
      i = l+1;
      j = ir;
      a = uar[l+1];
      b = sar[l+1];
      for (;;) {
	do i++; while (uar[i] < a);
	do j--; while (uar[j] > a);
	if (j < i) break;
	SWAP(uar[i], uar[j])
	  SWAP(sar[i], sar[j])
	  }
      uar[l+1] = uar[j];
      uar[j] = a;
      sar[l+1] = sar[j];
      sar[j] = b;
      jstack += 2;
      if (jstack >= n) exit_on_error("stack too small in sort");
      if (ir-i+1 >= j-l) {
	istack[jstack] = ir;
	istack[jstack-1] = i;
	ir = j-1;
      } else {
	istack[jstack] = j-1;
	istack[jstack-1] = l;
	l = i;
      }
    }
  }

}

void
sort_by_index (int n, double *uar, int *indx)
{
  int i, indxt, temp, ir=n-1, j, k, l=0, jstack=0, istack[n];
  double a;
  
  for (j = 0; j < n; j++) indx[j] = j;
  for (;;) {
    if (ir-l < 8) {
      for (j = l+1; j <= ir; j++) {
	indxt = indx[j];
	a = uar[indxt];
	for (i = j-1; i >= l; i--) {
	  if (uar[indx[i]] <= a) break;
	  indx[i+1] = indx[i];
	}
	indx[i+1] = indxt;
      }
      if (jstack == 0) break;
      ir = istack[jstack--];
      l  = istack[jstack--];
    } else {
      k = (l+ir) >> 1;
      SWAP(indx[k], indx[l+1]);
      if (uar[indx[l]] > uar[indx[ir]]) {
	SWAP(indx[l], indx[ir])
	  }
      if (uar[indx[l+1]] > uar[indx[ir]]) {
	SWAP(indx[l+1], indx[ir])
	  }
      if (uar[indx[l]] > uar[indx[l+1]]) {
	SWAP(indx[l], indx[l+1])
	  }
      i = l+1;
      j = ir;
      indxt = indx[l+1];
      a = uar[indxt];
      for (;;) {
	do i++; while (uar[indx[i]] < a);
	do j--; while (uar[indx[j]] > a);
	if (j < i) break;
	SWAP(indx[i], indx[j])
	  }
      indx[l+1] = indx[j];
      indx[j] = indxt;
      jstack += 2;
      if (jstack >= n) exit_on_error("stack too small in sort");
      if (ir-i+1 >= j-l) {
	istack[jstack] = ir;
	istack[jstack-1] = i;
	ir = j-1;
      } else {
	istack[jstack] = j-1;
	istack[jstack-1] = l;
	l = i;
      }
    }
  }

}

void
sort_by_indexi (int n, int *uar, int *indx)
{
  int i, indxt, temp, ir=n-1, j, k, l=0, jstack=0, istack[n];
  int a;
  
  for (j = 0; j < n; j++) indx[j] = j;
  for (;;) {
    if (ir-l < 2) {
      for (j = l+1; j <= ir; j++) {
	indxt = indx[j];
	a = uar[indxt];
	for (i = j-1; i >= l; i--) {
	  if (uar[indx[i]] <= a) break;
	  indx[i+1] = indx[i];
	}
	indx[i+1] = indxt;
      }
      if (jstack == 0) break;
      ir = istack[jstack--];
      l  = istack[jstack--];
    } else {
      k = (l+ir) >> 1;
      SWAP(indx[k], indx[l+1]);
      if (uar[indx[l]] > uar[indx[ir]]) {
	SWAP(indx[l], indx[ir])
	  }
      if (uar[indx[l+1]] > uar[indx[ir]]) {
	SWAP(indx[l+1], indx[ir])
	  }
      if (uar[indx[l]] > uar[indx[l+1]]) {
	SWAP(indx[l], indx[l+1])
	  }
      i = l+1;
      j = ir;
      indxt = indx[l+1];
      a = uar[indxt];
      for (;;) {
	do i++; while (uar[indx[i]] < a);
	do j--; while (uar[indx[j]] > a);
	if (j < i) break;
	SWAP(indx[i], indx[j])
	  }
      indx[l+1] = indx[j];
      indx[j] = indxt;
      jstack += 2;
      if (jstack >= n) exit_on_error("stack too small in sort");
      if (ir-i+1 >= j-l) {
	istack[jstack] = ir;
	istack[jstack-1] = i;
	ir = j-1;
      } else {
	istack[jstack] = j-1;
	istack[jstack-1] = l;
	l = i;
      }
    }
  }

}

void
sort_by_indexc (int n, char **uar, int *indx)
{
  int i, indxt, temp, ir=n-1, j, k, l=0, jstack=0, istack[n];
  char *a;
  
  for (j = 0; j < n; j++) indx[j] = j;
  for (;;) {
    if (ir-l < 8) {
      for (j = l+1; j <= ir; j++) {
	indxt = indx[j];
	a = ARRAY_ALLOC(strlen(uar[indxt]), char);
	a = strdup(uar[indxt]);
	for (i = j-1; i >= l; i--) {
	  if (strcmp(uar[indx[i]], a) <= 0) break;
	  indx[i+1] = indx[i];
	}
	free(a);
	indx[i+1] = indxt;
      }
      if (jstack == 0) break;
      ir = istack[jstack--];
      l  = istack[jstack--];
    } else {
      k = (l+ir) >> 1;
      SWAP(indx[k], indx[l+1]);
      if (strcmp(uar[indx[l]], uar[indx[ir]]) > 0) {
	SWAP(indx[l], indx[ir])
	  }
      if (strcmp(uar[indx[l+1]], uar[indx[ir]])> 0) {
	SWAP(indx[l+1], indx[ir])
	  }
      if (strcmp(uar[indx[l]], uar[indx[l+1]])> 0) {
	SWAP(indx[l], indx[l+1])
	  }
      i = l+1;
      j = ir;
      indxt = indx[l+1];
      a = ARRAY_ALLOC(strlen(uar[indxt]), char);
      a = strdup(uar[indxt]);
      for (;;) {
	do i++; while (strcmp(uar[indx[i]], a) < 0);
	do j--; while (strcmp(uar[indx[j]], a) > 0);
	if (j < i) break;
	SWAP(indx[i], indx[j])
	  }
      free(a);
      indx[l+1] = indx[j];
      indx[j] = indxt;
      jstack += 2;
      if (jstack >= n) exit_on_error("stack too small in sort");
      if (ir-i+1 >= j-l) {
	istack[jstack] = ir;
	istack[jstack-1] = i;
	ir = j-1;
      } else {
	istack[jstack] = j-1;
	istack[jstack-1] = l;
	l = i;
      }
    }
  }

}

void
reverse(int n, int *arr)
{

  int i, j, *rar;

  rar = ARRAY_ALLOC(n, int);

  j = n;
  for (i=0;i<n;i++) {
    j--;
    rar[j] = arr[i];
  }
  for (i=0;i<n;i++) {
    arr[i] = rar[i];
  }  

}

void
pearsn(double x[], double y[], int n, double *r, double *prob, double *z)
{
  double betai(double a, double b, double x);
  double erfcc(double x);
  int j;
  double yt,xt,t,df;
  double syy=0.0,sxy=0.0,sxx=0.0,ay=0.0,ax=0.0;
  
  for (j=0;j<n;j++) {
    ax += x[j];
    ay += y[j];
  }
  ax /= n;
  ay /= n;
  for (j=0;j<n;j++) {
    xt=x[j]-ax;
    yt=y[j]-ay;
    sxx += xt*xt;
    syy += yt*yt;
    sxy += xt*yt;
  }
  *r=sxy/(sqrt(sxx*syy)+TINY);
  *z=0.5*log((1.0+(*r)+TINY)/(1.0-(*r)+TINY));
  df=n-2;
  t=(*r)*sqrt(df/((1.0-(*r)+TINY)*(1.0+(*r)+TINY)));
  *prob=betai(0.5*df,0.5,df/(df+t*t));

}

void
spearman_corr(double data1[], double data2[], int n, int l, double *d, double *zd,
	      double *probd, double *rs, double *probrs)
{
  double betai(double a, double b, double x);
  void crank(int n, double w[], double *s);
  double erfcc(double x);
  void sort(int n, double *arr, double *brr);
  int j;
  double vard, t, sg, sf, fac, en3n, en, df, aved, *wksp1, *wksp2;

  wksp1 = ARRAY_ALLOC(n, double);
  wksp2 = ARRAY_ALLOC(n, double);
  for (j = 0; j < n; j++) {
    wksp1[j] = data1[j];
    wksp2[j] = data2[j];
  }

  sort(n, wksp1, wksp2);
  if (l == 1) {
    crankl(n, wksp1, &sf);
  } else {
    crank(n, wksp1, &sf);
  }
  sort(n, wksp2, wksp1);
  if (l == 1) {
    crankl(n, wksp2, &sg);
  } else {
    crank(n, wksp2, &sg);
  }
  *d = 0.0;
  for (j = 0; j < n; j++) {
    *d += SQR(wksp1[j] - wksp2[j]);
  }
  en = n;
  en3n = en*en*en-en;
  aved = en3n/6.0 - (sf+sg)/12.0;
  fac = (1.0-sf/en3n) * (1.0-sg/en3n);
  vard = ((en-1.0)*en*en*SQR(en+1.0)/36.0)*fac;
  *zd = (*d-aved)/sqrt(vard);
  *probd = erfcc(fabs(*zd)/1.4142136);
  *rs = (1.0-(6.0/en3n)*(*d+(sf+sg)/12.0))/sqrt(fac);
  fac = (*rs+1.0)*(1.0-(*rs));
  if (fac > 0.0) {
    t = (*rs)*sqrt((en-2.0)/fac);
    df = en-2.0;
    *probrs = betai(0.5*df,0.5,df/(df+t*t));
  } else
    *probrs = 0.0;
  free(wksp2);
  free(wksp1);

}

void
spearman_corr_nna(double data1tmp[], double data2tmp[], int ntmp, int l, double *d, double *zd,
		  double *probd, double *rs, double *probrs)
{
  int i, j;
  int n = 0;
  double *data1, *data2;
  double vard, t, sg, fac, en3n, en, df, aved;
  double betai(double a, double b, double x);
  void crank(int n, double w[], double *s);
  double erfcc(double x);
  void sort(int n, double *arr, double *brr);
    
  static double sf;
  static double *wksp1 = NULL;
  static double *wksp2 = NULL;
  static double *wksp3 = NULL;
  static double *wksp4 = NULL;
  
  // We drop the data with missing values
  data1 = ARRAY_ALLOC(ntmp, double);
  data2 = ARRAY_ALLOC(ntmp, double);
  for (i = 0; i < ntmp; i++) {
    double val = data1tmp[i] + data2tmp[i];
    if (! isnan(val)) {
      data1[n] = data1tmp[i];
      data2[n] = data2tmp[i];
      n++;
    }    
  }

  if (!wksp1) {
    wksp1 = ARRAY_ALLOC(n, double);
    wksp2 = ARRAY_ALLOC(n, double);
    wksp3 = ARRAY_ALLOC(n, double);
    wksp4 = ARRAY_ALLOC(n, double);
    for (j = 0; j < n; j++) {
      wksp2[j] = (double) j;
    }
    memcpy(wksp1, data1, sizeof(double)*n);
    sort(n, wksp1, wksp2);
    if (l == 1) {
      crankl(n, wksp1, &sf);
    } else {
      crank(n, wksp1, &sf);
    }
  }
  memcpy(wksp3, wksp1, sizeof(double)*n);

  for (j=0; j<n; j++) {
    wksp4[j] = data2[(int)(wksp2[j]+TINY)];
  }

  sort(n, wksp4, wksp3);
  if (l == 1) {
    crankl(n, wksp4, &sf);
  } else {
    crank(n, wksp4, &sf);
  }
  *d = 0.0;
  for (j = 0; j < n; j++) {
    *d += SQR(wksp3[j] - wksp4[j]);
  }
  en = n;
  en3n = en*en*en-en;
  aved = en3n/6.0 - (sf+sg)/12.0;
  fac = (1.0-sf/en3n) * (1.0-sg/en3n);
  vard = ((en-1.0)*en*en*SQR(en+1.0)/36.0)*fac;
  *zd = (*d-aved)/sqrt(vard);
  *probd = erfcc(fabs(*zd)/1.4142136);
  *rs = (1.0-(6.0/en3n)*(*d+(sf+sg)/12.0))/sqrt(fac);
  fac = (*rs+1.0)*(1.0-(*rs));
  if (fac > 0.0) {
    t = (*rs)*sqrt((en-2.0)/fac);
    df = en-2.0;
    *probrs = betai(0.5*df,0.5,df/(df+t*t));
  } else {
    *probrs = 0.0;
  }

  free(data1);
  free(data2);

}

double
betai(double a, double b, double x)
{
  double betacf(double a, double b, double x);
  double gammln(double xx);
  double bt, res;
  
  if (x < 0.0 || x > 1.0) {
    //exit_on_error("Bad x in routine betai");
    res = 10;
    return res;
  }
  if (x == 0.0 || x == 1.0) bt = 0.0;
  else
    bt = exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
  if (x < (a+1.0) / (a+b+2.0)) {
    res = bt*betacf(a,b,x)/a;
    if (isnan (res)) res=1.0;
    return res;
  } else {
    res = 1.0-bt*betacf(b,a,1.0-x)/b;
    if (isnan (res)) res=1.0;
    return res;
  }

}

double
betacf(double a, double b, double x)
{
  int m, m2;
  double aa, c, d, del, h, qab, qam, qap;
  
  qab = a + b;
  qap = a + 1.0;
  qam = a - 1.0;
  c = 1.0;
  d = 1.0-qab*x/qap;
  if (fabs(d) < FPMIN) d=FPMIN;
  d = 1.0/d;
  h = d;
  for (m = 1; m <= MAXIT; m++) {
    m2 = 2*m;
    aa = m*(b-m)*x/((qam+m2)*(a+m2));
    d = 1.0+aa*d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1.0+aa/c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d = 1.0+aa*d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1.0+aa/c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1.0/d;
    del = d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (m > MAXIT){
    //printf("a or b too big, or MAXIT too small in betacf\n");
    //printf("a = %f, b = %f\n", a, b);
  } 
  return h;

}

double
gammln(double xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  
  y = x = xx;
  tmp = x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser = 1.000000000190015;
  for (j = 0; j <= 5; j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);

}

void
rank(int n, int indx[], int irank[])
{

  int i;

  for (i=0;i<n;i++) irank[indx[i]]=i;

}


void
crank(int n, double w[], double *s)
{
  int j=0, ji, jt;
  double t, rank;
  
  *s = 0.0;
  while (j < n) {
    if (w[j+1] != w[j]) {
      w[j] = j;
      ++j;
    } else {
      for (jt = j+1; jt < n && w[jt] == w[j]; jt++);
      rank = 0.5*(j+jt-1);
      for (ji = j; ji <= (jt-1); ji++) w[ji] = rank;
      t = jt-j;
      *s += t*t*t-t;
      j = jt;
    }
  }

}

void
crankl(int n, double w[], double *s)
{
  int j=0, ji, jt;
  double t, rank;
  
  *s = 0.0;
  while (j < n) {
    if (w[j+1] != w[j]) {
      w[j] = j;
      ++j;
    } else {
      for (jt = j+1; jt < n && w[jt] == w[j]; jt++);
      rank = 0.5*(j+jt-1);
      for (ji = j; ji <= (jt-1); ji++) w[ji] = rank;
      t = jt-j;
      *s += log(t*t*t-t);
      j = jt;
    }
  }

}

double
erfcc(double x)
{
  double t,z,ans;
  
  z = fabs(x);
  t = 1.0/(1.0+0.5*z);
  ans = t*exp(-z*z-1.26551223+t*
	      (1.00002368+t*
	       (0.37409196+t*
		(0.09678418+t*
		 (-0.18628806+t*
		  (0.27886807+t*
		   (-1.13520398+t*
		    (1.48851587+t*
		     (-0.82215223+t*
		      0.17087277)))))))));
  return x >= 0.0 ? ans : 2.0-ans;

}

void
histogram(double *data, int n, int intervals, double bin, double min, int *hist)
{

  int i, b;
  double d;

  for (i = 0; i < intervals; i++) {
    hist[i] = 0; 
  }
  for (i = 0; i < n; i++) {
    b = 0;
    d = bin + min;
    while (data[i] > d) {
      d += bin;
      b++;
    }
    if (b >= intervals) {
      hist[intervals-1]++;
    } else {
      hist[b]++;
    }
  }

}

void
percentile(int *data, int n, int tot, double *prct)
{

  int i;
  double p = 0.0;

  for (i = 0; i < n; i++) {
    prct[i] = 0; 
  }
  for (i = 0; i < n; i++) {
    p += (double)data[i] / tot;
    prct[i] = p;
  }

}

void
shuffle(int *array, int n)
{
  int i, j, temp;
  
  for (i = n-1; i > 0; i--) {
    j = rand() % (i+1);
    SWAP (array[j], array[i]);
  }
  
}

int
factrl (int n)
{

  int fctrl=1;
  
  if (n < 0 ) exit_on_error("No factorial for negative numbers");
  if (n > 32) return exp(gammln(n+1.0));
  
  while (n > 1) {
    fctrl *= n;
    n--;
  }
  
  return fctrl;

}

int
permutation (int n, int k)
{

  int res=1;
  
  while (k--) res *= n--;
  
  return res;
  
}

int
choose (int n, int k)
{

  int res=1;
  int j=1;
  
  if (k > n || k < 0) return 0;
  if (n - k < k) k = n - k;
  
  while (j <= k) {
    res *= n--; 
    res /= j++;
  }
  
  return res;
    
}

void
avevar(double data[], int n, double *ave, double *var)
{
  int j;
  double s,ep;
  
  for (*ave=0.0,j=0;j<n;j++) *ave += data[j];
  *ave /= n;
  *var=ep=0.0;
  for (j=0;j<n;j++) {
    s=data[j]-(*ave);
    ep += s;
    *var += s*s;
  }
  *var=(*var-ep*ep/n)/(n-1);

}

void
bonferroni (double *data, double *cdata, int n)
{

  int i;

  for (i=0; i<n; i++) {
    cdata[i] = data[i]*n;
    if (cdata[i]>1) {
      cdata[i]=1.0;
    }
  }

}

void
fdr (int n, double *data, double *qval, double q,
	  double eta0, double *cutoff, int *nsig)
{

  int i, *o, *r, *s; 

  void sort_by_index(int m, double *arr, int *brr);
  void rank(int m, int indx[], int irank[]);

  s = ARRAY_ALLOC(n, int);
  o = ARRAY_ALLOC(n, int);
  r = ARRAY_ALLOC(n, int);
  sort_by_index(n,data,o);
  rank(n,o,r);

  for (i=0;i<n;i++) {
    qval[i]=eta0*n*data[i]/(r[i]+1);
  }
  qval[o[n-1]]=MIN(qval[o[n-1]],1);
  for (i=n-2;i>=0;i--) {
    qval[o[i]]=MIN(MIN(qval[o[i+1]],1),qval[o[i]]);
  }
  *nsig=0;
  for (i=0;i<n;i++) {
    if (qval[i]<=q) {
      *nsig += 1;
      s[i] = 1;
    } else {
      s[i] = 0;
    }
  } 
  *cutoff=0.0;
  if (*nsig!=0) {
    for (i=0;i<n;i++) {
      if (s[i]==1) {
	*cutoff=MAX(data[i],*cutoff);
      }
    }
  }
  for (i=0;i<n;i++) {
    if (data[i] <= *cutoff) {
      qval[i] = 1;
    } else {
      qval[i] = 0;
    }
  }

}

double
pythag(double a, double b)
{
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) {
    return absa*sqrt(1.0+SQR(absb/absa));
  } else {
    return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
  }
  
}

void
tqli(double d[], double e[], int n, double **z)
{
  double pythag(double a, double b);
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;

  for (i=1;i<n;i++) {
    e[i-1]=e[i];
  }
  e[n-1]=0.0;
  for (l=0;l<n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-1;m++) {
	dd=fabs(d[m])+fabs(d[m+1]);
	if ((fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
	if (iter++ == 30) exit_on_error("Too many iterations in tqli");
	g=(d[l+1]-d[l])/(2.0*e[l]);
	r=pythag(g,1.0);
	g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) {
	  f=s*e[i];
	  b=c*e[i];
	  e[i+1]=(r=pythag(f,g));
	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m]=0.0;
	    break;
	  }
	  s=f/r;
	  c=g/r;
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  d[i+1]=g+(p=s*r);
	  g=c*r-b;
	  for (k=0;k<n;k++) {
	    f=z[k][i+1];
	    z[k][i+1]=s*z[k][i]+c*f;
	    z[k][i]=c*z[k][i]-s*f;
	  }
	}
	if (r == 0.0 && i <= l) continue;
	d[l] -= p;
	e[l]=g;
	e[m]=0.0;
      }
    } while (m != l);
  }

}

void
tred2(double **a, int n, double d[], double e[])
{
  int l,k,j,i;
  double scale, hh,h,g,f;

  for (i=n-1;i>0;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 0) {
      for (k=0;k<=l;k++) {
	scale += fabs(a[i][k]);
      }
      if (scale == 0.0) {
	e[i]=a[i][l];
      } else {
	for (k=0;k<=l;k++) {
	  a[i][k] /= scale;
	  h += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
	e[i]=scale*g;
	h -= f*g;
	a[i][l]=f-g;
	f=0.0;
	for (j=0;j<=l;j++) {
	  a[j][i]=a[i][j]/h;
	  g=0.0;
	  for (k=0;k<=j;k++) {
	    g += a[j][k]*a[i][k];
	  }
	  for (k=j+1;k<=l;k++) {
	    g += a[k][j]*a[i][k];
	  }
	  e[j]=g/h;
	  f += e[j]*a[i][j];
	}
	hh=f/(h+h);
	for (j=0;j<=l;j++) {
	  f=a[i][j];
	  e[j]=g=e[j]-hh*f;
	  for (k=0;k<=j;k++) {
	    a[j][k] -= (f*e[k]+g*a[i][k]);
	  }
	}
      }
    } else {
      e[i]=a[i][l];
    }
    d[i]=h;
  }
  d[0]=0.0;
  e[0]=0.0;

  // Contents of this loop can be omitted if eigenvectors not
  // wanted except for statement d[i]=a[i][i];

  for (i=0;i<n;i++) {
    l=i-1;
    if (d[i]) {
      for (j=0;j<=l;j++) {
	g=0.0;
	for (k=0;k<=l;k++) {
	  g += a[i][k]*a[k][j];
	}
	for (k=0;k<=l;k++) {
	  a[k][j] -= g*a[k][i];
	}
      }
    }
    d[i]=a[i][i];
    a[i][i]=1.0;
    for (j=0;j<=l;j++) {
      a[j][i]=a[i][j]=0.0;
    }
  }

}

void
covsrt(double **covar, int ma, int ia[], int mfit)
{
    
  int i,j;
  double swap;
  
  for (j=0;j<ma-1;j++) {
    for (i=j+1;i<ma;i++) {
      covar[i][j]=0.0;
    }
  }
  for (i=0;i<mfit-1;i++) {
    for (j=i+1;j<mfit;j++) {
      if (ia[j]>ia[i]) {
	covar[ia[j]][ia[i]]=covar[i][j];
      }	else {
	covar[ia[i]][ia[j]]=covar[i][j];
      }
    }
  }
  swap=covar[0][0];
  for (j=0;j<ma;j++) {
    covar[0][j]=covar[j][j];
    covar[j][j]=0.0;
  }
  covar[ia[0]][ia[0]]=swap;
  for (j=1;j<mfit;j++) {
    covar[ia[j]][ia[j]]=covar[0][j];
  }
  for (j=1;j<ma;j++) {
    for (i=0;i<=j-1;i++) {
      covar[i][j]=covar[j][i];
    }
  }

}

double
probks(double alam)
{

  int j;
  double a2,fac=2.0,sum=0.0,term,termbf=0.0;
 
  a2 = -2.0*alam*alam;
  for (j=1;j<=100;j++) {
    term=fac*exp(a2*j*j);
    sum += term;
    if (fabs(term) <= EPSKS1*termbf || fabs(term) <= EPSKS2*sum) return sum;
    fac = -fac;
    termbf=fabs(term);
  }
  return 1.0;

}

void
zscore(double **matrix,
       int r,
       int c,
       int axis)
{

  int i, j, n;
  double mean, stdev;

  if (axis == 1) {

    for (i=0;i<r;i++) {
      mean = 0;
      stdev = 0;
      n = 0;
      for (j=0;j<c;j++) {
	if (! isnan(matrix[i][j])) {
	  mean += matrix[i][j];
	  stdev += matrix[i][j] * matrix[i][j];
	  n++;
	}
      }
      mean /= n;
      stdev = sqrt(stdev / n - mean * mean);
      for (j=0;j<c;j++) {
	if (! isnan(matrix[i][j])) {
	  matrix[i][j] = (matrix[i][j] - mean) / stdev;
	}
      }
    }  

  } else if (axis == 2) {

    for (i=0;i<c;i++) {
      mean = 0;
      stdev = 0;
      n = 0;
      for (j=0;j<r;j++) {
	if (! isnan(matrix[j][i])) {
	  mean += matrix[j][i];
	  stdev += matrix[j][i] * matrix[j][i];
	  n++;
	}
      }
      mean /= n;
      stdev = sqrt(stdev / n - mean * mean);
      for (j=0;j<r;j++) {
	if (! isnan(matrix[j][i])) {
	  matrix[j][i] = (matrix[j][i] - mean) / stdev;
	}
      }
    }  

  }

}
