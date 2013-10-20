/* randomd.c - double precision */
/* Routinen zur Erzeugung von Zufallszahlen */
/* Verwendet Routinen aus Numerical Recipes */
/* Volker Ahlers, Oktober 1998 */

#define _CRT_SECURE_NO_DEPRECATE
#include <math.h>
#include <time.h>

#include "randomd.h"
#include "nrutil.h"  /* Numerical Recipes */

#define RNG ran2

/* Zufallszahlen initialisieren */
void init_random(long *idum)
{
  time_t zeit;

  zeit = time((time_t *)NULL);
  *idum = -(long)(localtime(&zeit)->tm_sec);
}


/*---------------------------------------------------------------------*/  

/* Routinen aus Numerical Recipes, double-Versionen */
/* (double problematisch bei Zufallszahlen???) */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 2.3e-16
#define RNMX (1.0-EPS)


double ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 2.3e-16
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX



double gasdev(long *idum)
{
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;

  if  (iset == 0) {
    do {
      v1=2.0*RNG(idum)-1.0;
      v2=2.0*RNG(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}


double expdev(long *idum)
{
  double dum;

  do
    dum=RNG(idum);
  while (dum == 0.0);
  return -log(dum);
}

/*----------------------------------------------------------------*/

/* Gauﬂverteilte Zufallszahl, E(x)=mean, sqrt(Var(x))=stddev */
double gaussrandom(double mean,double stddev,long *idum)
{
  if(stddev <= 0)
    nrerror("gaussrandom: stddev must be positive");
  return stddev*gasdev(idum)+mean;
}


/* Im Intervall (left,right) gleichverteilte Zufallszahl */
double unirandom(double left,double right,long *idum)
{
  if(left >= right)
    nrerror("idrandom: left must be less than right");
  return RNG(idum)*(right-left) + left;
}


/* Exponentiell verteilte Zufallszahl, p(x)=lambda*exp(-lambda*x) */
/* E(x)=1/lambda */
/*double exprandom(double lambda,long *idum)*/
/*{*/
  /*if(lambda <= 0)*/
    /*nrerror("exprandom: lambda must be positive");*/
  /*return expdev(idum)/lambda;*/
/*}*/
