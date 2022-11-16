/*********************************************************************
 Copyright 2004 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: solution.c,v 1.21 2009/04/02 21:41:11 neuhausi Exp $
**********************************************************************/

#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <unistd.h>
#include <math.h>
#define EXTERN
#include "solution.h"

void
ludcmp(double **a, int n, int *indx, double *d)
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double vv[n];

  *d=1.0;
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) exit_on_error("Singular matrix in routine ludcmp");
    vv[i]=1.0/big;
  }
  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
      sum=a[i][j];
      for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<n;i++) {
      sum=a[i][j];
      for (k=0;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<n;i++) a[i][j] *= dum;
    }
  }

}

void
lubksb(double **a, int n, int *indx, double *b)
{
  int i,ii=0,ip,j;
  double sum;

  for (i=0;i<n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    //if (ii)
      for (j=ii;j<i;j++) sum -= a[i][j]*b[j];
      //else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n-1;i>=0;i--) {
    sum=b[i];
    for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }

}

void
svdcmp(double **a, int m, int n, double *w, double **v)
{

  int flag,i,its,j,jj,k,l,nm;
  double anorm,c,f,g,h,s,scale,x,y,z,rv1[n];

  g=scale=anorm=0.0;
  for (i=0;i<n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;

    if (i < m) {
      for (k=i;k<m;k++) {
	scale += fabs(a[k][i]);
      }
      if (scale) {
	for (k=i;k<m;k++) {
	  a[k][i] /= scale;
	  s += a[k][i]*a[k][i];
	}
	f=a[i][i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][i]=f-g;
	for (j=l;j<n;j++) {
	  for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
	  f=s/h;
	  for (k=i;k<m;k++) a[k][j] += f*a[k][i];
	}
	for (k=i;k<m;k++) {
	  a[k][i] *= scale;
	}
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i < m && i != n-1) {
      for (k=l;k<n;k++) {
	scale += fabs(a[i][k]);
      }
      if (scale) {
	for (k=l;k<n;k++) {
	  a[i][k] /= scale;
	  s += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][l]=f-g;
	for (k=l;k<n;k++) {
	  rv1[k]=a[i][k]/h;
	}
	for (j=l;j<m;j++) {
	  for (s=0.0,k=l;k<n;k++) {
	    s += a[j][k]*a[i][k];
	  }
	  for (k=l;k<n;k++) {
	    a[j][k] += s*rv1[k];
	  }
	}
	for (k=l;k<n;k++) {
	  a[i][k] *= scale;
	}
      }
    }
    anorm=MAX(anorm,((double)fabs(w[i])+(double)fabs(rv1[i])));
  }

  for (i=n-1;i>=0;i--) {
    if (i < n-1) {
      if (g) {
	for (j=l;j<n;j++) {
	  v[j][i]=(a[i][j]/a[i][l])/g;
	}
	for (j=l;j<n;j++) {
	  for (s=0.0,k=l;k<n;k++) {
	    s += a[i][k]*v[k][j];
	  }
	  for (k=l;k<n;k++) {
	    v[k][j] += s*v[k][i];
	  }
	}
      }
      for (j=l;j<n;j++) {
	v[i][j]=v[j][i]=0.0;
      }
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }

  for (i=MIN(m,n)-1;i>=0;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<n;j++) {
      a[i][j]=0.0;
    }
    if (g) {
      g=1.0/g;
      for (j=l;j<n;j++) {
	for (s=0.0,k=l;k<m;k++) {
	  s += a[k][i]*a[k][j];
	}
	f=(s/a[i][i])*g;
	for (k=i;k<m;k++) {
	  a[k][j] += f*a[k][i];
	}
      }
      for (j=i;j<m;j++) {
	a[j][i] *= g;
      }
    } else {
      for (j=i;j<m;j++) {
	a[j][i]=0.0;
      }
    }
    ++a[i][i];
  }

  for (k=n-1;k>=0;k--) {
    for (its=0;its<30;its++) {
      flag=1;
      for (l=k;l>=0;l--) {
	nm=l-1;
	if ((double)(fabs(rv1[l])+anorm) == anorm) {
	  flag=0;
	  break;
	}
	if ((double)(fabs(w[nm])+anorm) == anorm) {
	  break;
	}
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  if ((double)(fabs(f)+anorm) != anorm) break;
	  g=w[i];
	  h=pythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=0;j<m;j++) {
	    y=a[j][nm];
	    z=a[j][i];
	    a[j][nm]=y*c+z*s;
	    a[j][i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=0;j<n;j++) {
	    v[j][k] = -v[j][k];
	  }
	}
	break;
      }

      if (its == 29) {
	exit_on_error("no convergence in 30 svdcmp iterations");
      }

      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g = g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=0;jj<n;jj++) {
	  x=v[jj][j];
	  z=v[jj][i];
	  v[jj][j]=x*c+z*s;
	  v[jj][i]=z*c-x*s;
	}
	z=pythag(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=0;jj<m;jj++) {
	  y=a[jj][j];
	  z=a[jj][i];
	  a[jj][j]=y*c+z*s;
	  a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }

  //if (debug) {
  //  printf (" w\n");
  //  print_array_d(w, n);	
  //  printf (" v\n");
  //  for (i=0;i<n;i++) {
  //    print_array_d(v[i], n);	
  //  }
  //  printf (" a\n");
  //  for (i=0;i<m;i++) {
  //    print_array_d(a[i], n);	
  //  }
  //}

}

void
svbksb(double **u, double *w, double **v, int m, int n, double *b, double *x)
{

  int jj,j,i;
  double s,tmp[n];

  for (j=0;j<n;j++) {
    s=0.0;
    if (w[j]) {
      for (i=0;i<m;i++) s += u[i][j]*b[i];
      s /= w[j];
    }
    tmp[j]=s;
  }
  for (j=0;j<n;j++) {
    s=0.0;
    for (jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
    x[j]=s;
  }

}

void
qrdcmp(double **a, int n, double *c, double *d, int *sing)
{
  int i,j,k;
  double scale,sigma,sum,tau;

  *sing=0;
  for (k=0;k<n-1;k++) {
    scale=0.0;
    for (i=k;i<n;i++) scale=MAX(scale,fabs(a[i][k]));
    if (scale == 0.0) {
      *sing=1;
      c[k]=d[k]=0.0;
    } else {
      for (i=k;i<n;i++) a[i][k] /= scale;
      for (sum=0.0,i=k;i<n;i++) sum += SQR(a[i][k]);
      sigma=SIGN(sqrt(sum),a[k][k]);
      a[k][k] += sigma;
      c[k]=sigma*a[k][k];
      d[k] = -scale*sigma;
      for (j=k+1;j<n;j++) {
	for (sum=0.0,i=k;i<n;i++) sum += a[i][k]*a[i][j];
	tau=sum/c[k];
	for (i=k;i<n;i++) a[i][j] -= tau*a[i][k];
      }
    }
  }
  d[n]=a[n][n];
  if (d[n] == 0.0) *sing=1;

}

void
qrsolv(double **a, int n, double *c, double *d, double *b)
{

  void rsolv(double **a, int n, double *d, double *b);
  int i,j;
  double sum,tau;

  for (j=0;j<n-1;j++) {
    for (sum=0.0,i=j;i<n;i++) sum += a[i][j]*b[i];
    tau=sum/c[j];
    for (i=j;i<n;i++) b[i] -= tau*a[i][j];
  }
  rsolv(a,n,d,b);

}

void
rsolv(double **a, int n, double *d, double *b)
{
  int i,j;
  double sum;

  b[n] /= d[n];
  for (i=n-2;i>=0;i--) {
    for (sum=0.0,j=i+1;j<n;j++) sum += a[i][j]*b[j];
    b[i]=(b[i]-sum)/d[i];
  }

}

void
qrupdt(double **r, double **qt, int n, double *u, double *v)
{
  void rotate(double **r, double **qt, int n, int i, double a, double b);
  int i,j,k;

  for (k=n-1;k>=0;k--) {
    if (u[k]) break;
  }
  if (k < 1) k=1;
  for (i=k-2;i>=0;i--) {
    rotate(r,qt,n,i,u[i],-u[i+1]);
    if (u[i] == 0.0) u[i]=fabs(u[i+1]);
    else if (fabs(u[i]) > fabs(u[i+1]))
      u[i]=fabs(u[i])*sqrt(1.0+SQR(u[i+1]/u[i]));
    else u[i]=fabs(u[i+1])*sqrt(1.0+SQR(u[i]/u[i+1]));
  }
  for (j=0;j<n;j++) r[1][j] += u[1]*v[j];
  for (i=0;i<k;i++)
    rotate(r,qt,n,i,r[i][i],-r[i+1][i]);

}

void
rotate(double **r, double **qt, int n, int i, double a, double b)
{

  int j;
  double c,fact,s,w,y;

  if (a == 0.0) {
    c=0.0;
    s=(b >= 0.0 ? 1.0 : -1.0);
  } else if (fabs(a) > fabs(b)) {
    fact=b/a;
    c=SIGN(1.0/sqrt(1.0+(fact*fact)),a);
    s=fact*c;
  } else {
    fact=a/b;
    s=SIGN(1.0/sqrt(1.0+(fact*fact)),b);
    c=fact*s;
  }
  for (j=i;j<n;j++) {
    y=r[i][j];
    w=r[i+1][j];
    r[i][j]=c*y-s*w;
    r[i+1][j]=s*y+c*w;
  }
  for (j=0;j<n;j++) {
    y=qt[i][j];
    w=qt[i+1][j];
    qt[i][j]=c*y-s*w;
    qt[i+1][j]=s*y+c*w;
  }

}

void
cholsl(double **a, int n, double *p, double *b, double *x)
{

  int i,k;
  double sum;

  for (i=0;i<n;i++) {
    for (sum=b[i],k=i-2;k>=0;k--) sum -= a[i][k]*x[k];
    x[i]=sum/p[i];
  }
  for (i=n-1;i>=0;i--) {
    for (sum=x[i],k=i+1;k<n;k++) sum -= a[k][i]*x[k];
    x[i]=sum/p[i];
  }

}

void
choldc(double **a, int n, float *p)
{

  int i,j,k;
  double sum;

  for (i=0;i<n;i++) {
    for (j=i;j<n;j++) {
      sum=a[i][j];
      for (k=i-1;k>=0;k--) {
	sum -= a[i][k]*a[j][k];
      }
      if (i == j) {
	if (sum <= 0.0)
	  exit_on_error("choldc failed");
	p[i]=sqrt(sum);
      } else a[j][i]=sum/p[i];
    }
  }

}

void
cholroot(double **a, int n, float *p)
{

  int i, ii, j, jj, zeros[n], cr=0;
  double chol[n][n], **cchol;
  double zero[n];

  cchol = ARRAY_ALLOC(n, double *);
  for(i=0;i<n;i++) {
    cchol[i] = ARRAY_ALLOC(n, double);
    for(j=0;j<n;j++) {
      cchol[i][j] = 0;
    }
  }

  // Initialize zero array
  for(i=0;i<n;i++) {
    zero[i]=0;
  }

  // Sum all the columns
  for(i=0;i<n;i++) {
    for(j=0;j<n;j++) {
      zero[j] += fabs(a[i][j]);
    }
  }

  // Count the number of non zero rows/columns
  for(i=0;i<n;i++) {
    if (zero[i]>EPS) {
      cr++;
    }
  }

  // Copy non zero rows/columns toa new matrix
  ii = 0;
  for(i=0;i<n;i++) {
    jj=0;
    for(j=0;j<n;j++) {
      if (zero[j]>0) {
	zeros[j] = jj;
	cchol[ii][jj] = a[i][j];
	jj++;
      } else {
	zeros[j] = -1;
      }
    }
    if (zero[i]>0) {
      ii++;
    }
  }

  // Calculate cholesky
  choldc(cchol,cr,p);

  // Copy back the data into the a matrix
  for(i=0;i<n;i++) {
    if (zeros[i] > -1) {
      for(j=0;j<n;j++) {
	if (zeros[j] > -1) {
	  a[i][j] = cchol[zeros[i]][zeros[j]];
	}
      }
    }
  }

  for (i=0;i<n;i++)
     for (j=0;j<n;j++) chol[i][j]=a[i][j];

  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      a[i][j]=((i > j) ? chol[i][j] : (i == j ? p[zeros[i]] : 0.0));
      if (i > j) chol[i][j]=a[i][j];
      else chol[i][j]=(i == j ? p[zeros[i]] : 0.0);
    }
  }

  free(cchol);

}

int
gaussj(double **a, int n, double **b, int m)
{
  int indxc[n],indxr[n],ipiv[n];
  int i,icol,irow,j,k,l,ll;
  double big,dum,pivinv,temp;

  for (j=0;j<n;j++) {
    ipiv[j]=0;
  }

  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++) {
      if (ipiv[j] != 1) {
	for (k=0;k<n;k++) {
	  if (ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big=fabs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  }
	}
      }
    }
    ++(ipiv[icol]);
    if (irow != icol) {
      // return 0;
      // The following lines are returning an error
      for (l=0;l<n;l++) {
	SWAP(a[irow][l],a[icol][l]);
      }
      for (l=0;l<m;l++) {
	SWAP(b[irow][l],b[icol][l]);
      }
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) {
      return 0;
    }
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=0;l<n;l++) {
      a[icol][l] *= pivinv;
    }
    for (l=0;l<m;l++) {
      b[icol][l] *= pivinv;
    }
    for (ll=0;ll<n;ll++) {
      if (ll != icol) {
	dum=a[ll][icol];
	a[ll][icol]=0.0;
	for (l=0;l<n;l++) {
	  a[ll][l] -= a[icol][l]*dum;
	}
	for (l=0;l<m;l++) {
	  b[ll][l] -= b[icol][l]*dum;
	}
      }
    }
  }
  for (l=n-1;l>=0;l--) {
    if (indxr[l] != indxc[l]) {
      for (k=0;k<n;k++) {
	SWAP(a[k][indxr[l]],a[k][indxc[l]]);
      }
    }
  }

  return 1;

}

int
mrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[],
       int ma, double **covar, double **alpha, double *chisq,
       void (*funcs)(double, double [], double *, double [], int), double *alamda)
{

  void covsrt(double **covar, int ma, int ia[], int mfit);
  int gaussj(double **a, int n, double **b, int m);
  void mrqcof(double x[], double y[], double sig[], int ndata, double a[],
	      int ia[], int ma, double **alpha, double beta[], double **chisq,
	      void (*funcs)(double, double [], double *, double [], int));
  int j,k,l;
  static int mfit;
  static double ochisq,*atry,*beta,*da,**oneda;

  if (*alamda < 0.0) {
    atry=ARRAY_ALLOC(ma, double);
    beta=ARRAY_ALLOC(ma, double);
    da=ARRAY_ALLOC(ma, double);
    for (mfit=0,j=0;j<ma;j++) {
      if (ia[j]) {
	mfit++;
      }
    }
    oneda=ARRAY_ALLOC(mfit, double *);
    for (j=0;j<mfit;j++) {
      oneda[j]=ARRAY_ALLOC(1, double);
    }
    *alamda=0.001;
    mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,&chisq,funcs);
    ochisq=(*chisq);
    for (j=0;j<ma;j++) {
      atry[j]=a[j];
    }
  }
  for (j=0;j<mfit;j++) {
    for (k=0;k<mfit;k++) {
      covar[j][k]=alpha[j][k];
    }
    covar[j][j]=alpha[j][j]*(1.0+(*alamda));
    oneda[j][0]=beta[j];
  }
  if (! gaussj(covar,mfit,oneda,1)) {
    return 0;
  }
  for (j=0;j<mfit;j++) {
    da[j]=oneda[j][0];
  }
  if (*alamda == 0.0) {
    covsrt(covar,ma,ia,mfit);
    covsrt(alpha,ma,ia,mfit);
    free(oneda);
    free(da);
    free(beta);
    free(atry);
    return 1;
  }
  for (j=0,l=0;l<ma;l++) {
    if (ia[l]) {
      atry[l]=a[l]+da[j++];
    }
  }
  mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,&chisq,funcs);
  if (*chisq < ochisq) {
    *alamda *= 0.1;
    ochisq=(*chisq);
    for (j=0;j<mfit;j++) {
      for (k=0;k<mfit;k++) {
	alpha[j][k]=covar[j][k];
      }
      beta[j]=da[j];
    }
    for (l=0;l<ma;l++) {
      a[l]=atry[l];
    }
  } else {
    *alamda *= 10.0;
    *chisq=ochisq;
  }

  return 1;

}

void
mrqcof(double x[], double y[], double sig[], int ndata, double a[], int ia[],
       int ma, double **alpha, double beta[], double **chisq,
       void (*funcs)(double, double [], double *, double [], int))
{
	
  int i,j,k,l,m,mfit=0;
  double ymod,wt,sig2i,dy,dyda[ma];

  for (j=0;j<ma;j++) {
    if (ia[j]) {
      mfit++;
    }
  }
  for (j=0;j<mfit;j++) {
    for (k=0;k<=j;k++) {
      alpha[j][k]=0.0;
    }
    beta[j]=0.0;
  }
  **chisq=0.0;
  for (i=0;i<ndata;i++) {
    (*funcs)(x[i],a,&ymod,dyda,ma);
    sig2i=1.0/(sig[i]*sig[i]);
    dy=y[i]-ymod;
    for (j=-1,l=0;l<ma;l++) {
      if (ia[l]) {
	wt=dyda[l]*sig2i;
	for (++j,k=0,m=0;m<=l;m++) {
	  if (ia[m]) {
	    alpha[j][k++] += wt*dyda[m];
	  }
	}
	beta[j] += dy*wt;
      }
    }
    **chisq += dy*dy*sig2i;
  }
  for (j=1;j<mfit;j++) {
    for (k=0;k<j;k++) {
      alpha[k][j]=alpha[j][k];
    }
  }

}
