/* vecmatd.c - Routinen zur Vektor- und Matrixrechnung */
/*             double precision                        */
/* Arbeiten mit Arraydarstellung aus Numerical Recipes */
/* Volker Ahlers, November 1998                        */
/* cancelled mult_mat_mat because of compiling problem, Rüdiger 2004 */

#include <stdio.h>
#include <math.h>

#include "vecmatd.h"


/* Vektor drucken */
/* x besitzt dim Zeilen */
void print_vec(double *x,int dim)
{
  int i;

  for(i=1;i<=dim;i++)
    printf("%8.3g\n",x[i]);
}


/* Matrix drucken */
/* A besitzt dim1 Zeilen und dim2 Spalten */
void print_mat(double **a,int dim1,int dim2)
{
  int i,j;

  for(i=1;i<=dim1;i++) {
    for(j=1;j<=dim2;j++)
      printf(" %8.3g",a[i][j]);
    printf("\n");
  }
}


/* Einheitsmatrix erzeugen: A = I */
/* A besitzt dim Zeilen/Spalten */
void set_mat_id(double **a,int dim)
{
  int i,j;
  
  for(i=1;i<=dim;i++) {
    for(j=1;j<i;j++)
      a[i][j] = 0.;
    a[i][i] = 1.;
    for(j=i+1;j<=dim;j++)
      a[i][j] = 0.;
  }
}


/* Vektor kopieren: y = x */
/* x besitzt dim Zeilen */
void copy_vec(double *y,double *x,int dim)
{
  int i;

  for(i=1;i<=dim;i++)
    y[i] = x[i];
}


/* Matrix kopieren: B = A */
/* A besitzt dim1 Zeilen und dim2 Spalten */
void copy_mat(double **b,double **a,int dim1,int dim2)
{
  int i,j;

  for(i=1;i<=dim1;i++)
    for(j=1;j<=dim2;j++)
      b[i][j] = a[i][j];
}


/* Skalarprodukt: a = <x,y> */
/* x besitzt dim Zeilen, Funktionswert ist a */
double skalprod(double *x,double *y,int dim)
{
  int i;
  double prod=0.;

  for(i=1;i<=dim;i++)
    prod += x[i]*y[i];
  return prod;
}


/* Norm eines Vektors: |x| */
/* x besitzt dim Zeilen, Funktionswert ist |x| */
double norm_vec(double *x,int dim)
{
  int i;
  double norm;

  norm = x[1]*x[1];
  for(i=2;i<=dim;i++)
    norm += x[i]*x[i];
  return sqrt(norm);
}


/* Vektor renormieren: x -> x/|x| */
/* x besitzt dim Zeilen, Funktionswert ist |x| */
double renorm_vec(double *x,int dim)
{
  int i;
  double norm;

  norm = x[1]*x[1];
  for(i=2;i<=dim;i++)
    norm += x[i]*x[i];
  norm = sqrt(norm);
  for(i=1;i<=dim;i++)
    x[i] /= norm;
  return norm;
}


/* Matrix-Vektor-Multiplikation: y = Ax */
/* A besitzt dim1 Zeilen und dim2 Spalten */
void mult_mat_vec(double *y,double **a,double *x,int dim1,int dim2)
{
  int i,j;

  for(i=1;i<=dim1;i++) {
    y[i] = a[i][1]*x[1];
    for(j=2;j<=dim2;j++) {
	y[i] += a[i][j]*x[j];
    }
  }
}


/* Matrix-Matrix-Multiplikation: C = AB */
/* A besitzt dim1 Zeilen und dim2 Spalten, B dim2 Zeilen und dim3 Spalten */
/* => C besitzt dim1 Zeilen und dim3 Spalten */
/*void mult_mat_mat(double **c,double **a,double **b,int dim1,int dim2,int
dim3)*/
/*{*/
  /*int i,j,k;*/

  /*for(i=1;i<=dim1;i++)*/
    /*for(j=1;j<=dim3;j++) {*/
      /*c[i][j] = a[i][1]*b[1][j];*/
      /*for(k=2;k<=dim2;k++)*/
	/*c[i][j] += a[i][k]*b[k][j];*/
    /*}*/
/*}*/


/* ALT: Matrix-Matrix-Multiplikation: C = AB */
/* A besitzt dim1 Zeilen und dim2 Spalten, B dim2 Zeilen und Spalten */
/* => C besitzt dim1 Zeilen und dim2 Spalten */
void mult_mat_mat_alt(double **c,double **a,double **b,int dim1,int dim2)
{
  int i,j,k;

  for(i=1;i<=dim1;i++)
    for(j=1;j<=dim2;j++) {
      c[i][j] = a[i][1]*b[1][j];
      for(k=2;k<=dim2;k++)
	c[i][j] += a[i][k]*b[k][j];
    }
}


/* klass. Gram-Schmidt-Verfahren fuer m*n-Matrix Q       */
/* Diagonalelemente von R werden in r[1...n] gespeichert */
/* Hilfsvektor qtemp[1...m] wird benötigt                */
void cgs(double **q,double *qtemp,double *r,int m,int n)
{
  int i,k,l;
  double prod;

  for(k=1;k<=n;k++) {
    for(l=1;l<=m;l++)
      qtemp[l] = q[k][l];
    for(i=1;i<k;i++) {
      prod = q[i][1]*qtemp[1];
      for(l=2;l<=m;l++)
	prod += q[i][l]*qtemp[l];
      for(l=1;l<=m;l++)
	q[k][l] -= prod*q[i][l];
    }
    prod = q[k][1]*q[k][1];
    for(l=2;l<=m;l++)
      prod += q[k][l]*q[k][l];
    prod = sqrt(prod);
    r[k] = prod;
    for(l=1;l<=m;l++)
      q[k][l] /= prod;
  }
}


/* modifiziertes Gram-Schmidt-Verfahren fuer m*n-Matrix Q */
/* Diagonalelemente von R werden in r[1...n] gespeichert  */
/* zeilenorientiert                                       */
void mgsr(double **q,double *r,int m,int n)
{
  int j,k,l;
  double prod;

  for(k=1;k<=n;k++) {
    prod = q[k][1]*q[k][1];
    for(l=2;l<=m;l++)
      prod += q[k][l]*q[k][l];
    prod = sqrt(prod);
    r[k] = prod;
    for(l=1;l<=m;l++)
      q[k][l] /= prod;
    for(j=k+1;j<=n;j++) {
      prod = q[k][1]*q[j][1];
      for(l=2;l<=m;l++)
	prod += q[k][l]*q[j][l];
      for(l=1;l<=m;l++)
	q[j][l] -= prod*q[k][l];
    }
  }
}


/* modifiziertes Gram-Schmidt-Verfahren fuer m*n-Matrix Q */
/* ohne Speichern der Diagonalelemente */
/* zeilenorientiert */
void mgsr_solo(double **q,int m,int n)
{
  int j,k,l;
  double prod;

  for(k=1;k<=n;k++) {
    prod = q[k][1]*q[k][1];
    for(l=2;l<=m;l++)
      prod += q[k][l]*q[k][l];
    prod = sqrt(prod);
    for(l=1;l<=m;l++)
      q[k][l] /= prod;
    for(j=k+1;j<=n;j++) {
      prod = q[k][1]*q[j][1];
      for(l=2;l<=m;l++)
	prod += q[k][l]*q[j][l];
      for(l=1;l<=m;l++)
	q[j][l] -= prod*q[k][l];
    }
  }
}


/* modifiziertes Gram-Schmidt-Verfahren fuer m*n-Matrix Q */
/* Diagonalelemente von R werden in r[1...n] gespeichert  */
/* zeilenorientiert, Summennorm                           */
void mgsr_sum(double **q,double *r,int m,int n)
{
  int j,k,l;
  double prod;

  for(k=1;k<=n;k++) {
    prod = fabs(q[k][1]);
    for(l=2;l<=m;l++)
      prod += fabs(q[k][l]);
    r[k] = prod;
    for(l=1;l<=m;l++)
      q[k][l] /= prod;
    for(j=k+1;j<=n;j++) {
      prod = q[k][1]*q[j][1];
      for(l=2;l<=m;l++)
	prod += q[k][l]*q[j][l];
      for(l=1;l<=m;l++)
	q[j][l] -= prod*q[k][l];
    }
  }
}


/* Orthonormale Matrix A ohne Nullkomponenten erzeugen */
/* A besitzt dim Zeilen/Spalten */
void set_mat_orthonorm(double **a,int dim)
{
  int i,j;

  for(j=1;j<=dim;j++)
    a[1][j] = 1.;
  for(i=2;i<=dim;i++) {
    for(j=1;j<i;j++)
      a[i][j] = 1.;
    a[i][i] = -1.;
    for(j=i+1;j<=dim;j++)
      a[i][j] = 1.;
  }
  mgsr_solo(a,dim,dim);
}
