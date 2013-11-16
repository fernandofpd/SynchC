/* swarmsyncRefRan.c - double precision */
/* unit free */
/* random re-orientation at boundary */
/* Verwendet Routinen aus Numerical Recipes */
/* Ruediger Zillmer, May 2013 */
/* Fernando Perez Diaz, November 2013 */

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "nrutil.h"
#include "vecmatd.h"
#include "randomd.h"
#include "options.h"

/*#define DIM 4*/
/*#define NLYAP 2*/
#define Pi         2.*asin(1.)

/* ---------- parameter for random number generator -------------------- */

typedef unsigned long uint32;

#define NORM 4294967295
#define N (624)
#define M (397)
#define K (0x9908B0DFU)
#define hiBit(u) ((u) & 0x80000000U)
#define loBit(u) ((u) & 0x00000001U)
#define loBits(u) ((u) & 0x7FFFFFFFU)
#define mixBits(u,v) (hiBit(u)|loBits(v))

static uint32 state[N+1];
static uint32 *next;
static int left = -1;

/* ---------------------------------------------------- */

void maxfind(int *p,double *Gam);
void sort(unsigned long n, double a[], int b[]);
/* random number generator */
void seedMT(uint32 seed);
double ranMT(void);
/* find nearest neighbor of unit k */
int nearest(double *px,double *py,int k);
double incremod(double r,double v,double D);
double orderpar(double *Gam);
void bhit(double *x,double *y,double *vx,double *vy,double *Gam,double *t);
void printOrderpar(FILE *output, double *Gam);

/* Parameter: */
double tau=1.;
double dphi=0.1;
/* epsilon, unit number S, box size L, velocity Vel */
double eps=0.1;
int S=20;
double L=400;
double Vel=0.1;

int main(int argc,char **argv)
{
    int i,j,k,seed1=7,n=200,p,nn,fin=1;
    uint32 seed=5;
    double t,time=0,phase,zdum;
    /* neuron phase Gam, motion angle Phi */ 
    /* xy-position Pxy, xy-velocity Vxy */
    double *Gam,*Phi,*Px,*Py,*Vx,*Vy;
    int *P;    
    FILE *out1;

    if(opt_flag(&argc,argv,"-h")) {
        printf("\nOptions: -n\n\n");
        exit(0);
    }

    opt_int(&argc,argv,"-n",&n);
    opt_int(&argc,argv,"-N",&S);
    opt_int(&argc,argv,"-s",&seed1);
    opt_double(&argc,argv,"-dp",&dphi);
    opt_double(&argc,argv,"-V",&Vel);
    opt_double(&argc,argv,"-e",&eps);
    opt_double(&argc,argv,"-L",&L);

    out1=fopen("dat","w"); fclose(out1); out1=fopen("dat","a");
    
    P = ivector(1,S);
    Gam = dvector(1,S);
    Phi = dvector(1,S);
    Px = dvector(1,S);
    Py = dvector(1,S);
    Vx = dvector(1,S);
    Vy = dvector(1,S);

    /* initialize random numbers */
    seed=(uint32)(seed1); seedMT(seed);
    /* random initial values */
    for(i=1;i<=S;i++) {
        Gam[i]=ranMT(); Phi[i]=2*Pi*ranMT();
        Px[i]=L*ranMT(); Py[i]=L*ranMT();
        Vx[i]=Vel*cos(Phi[i]); Vy[i]=Vel*sin(Phi[i]);
    } 

    for(i=1;i<=n;i++) {
        fin=1;
        /* index of next firing unit p, time t */
        maxfind(&p,Gam); t = 1-Gam[p];
        /* update position and velocities */        
        bhit(Px,Py,Vx,Vy,Gam,&t);

        /* when reference unit 1 fires print out order parameter */
        if(p==1) printOrderpar(out1, Gam);
        Gam[p]=0;
        /* the loop for the case that firing triggers other firings */
        while(fin > 0) {
            nn=nearest(Px,Py,p);                
            Gam[nn] *= (1+eps);                
            if(Gam[nn] >= 1) {
                p = nn; Gam[nn]=0;
                if(p==1) printOrderpar(out1, Gam);
            }
            else fin=0;
        }        
    } 

    fclose(out1);

    free_ivector(P,1,S);
    free_dvector(Gam,1,S);

    return 0;
}

void maxfind(int *p,double *Gam)
{
  int i;
  
  *p=1; 
  for(i=2;i<=S;i++) if(Gam[i] > Gam[*p]) *p=i; 
}

void sort(unsigned long n, double a[], int b[])
{
	unsigned long i,j,inc;
        int w;
	double v;
	inc=1;
	do {
		inc *= 3;
		inc++;
	} while (inc <= n);
	do {
		inc /= 3;
		for (i=inc+1;i<=n;i++) {
			v=a[i]; w=b[i];
			j=i;
			while (a[j-inc] > v) {
				a[j]=a[j-inc]; b[j]=b[j-inc];
				j -= inc;
				if (j <= inc) break;
			}
			a[j]=v; b[j]=w;
		}
	} while (inc > 1);
}

/* random number generator cokus.c */

void seedMT(uint32 seed)
{
  register uint32 x = (seed | 1U) & 0xFFFFFFFFU,*s=state;
  register int j;
  
  for(left=0,*s++=x,j=N;--j;*s++ = (x*=69069U) & 0xFFFFFFFFU);
}

uint32 reloadMT(void)
{
  register uint32 *p0=state, *p2=state+2, *pM=state+M, s0,s1;
  register int j;
  
  if(left < -1) seedMT(4357U);
  left=N-1; next=state+1;
  for(s0=state[0],s1=state[1],j=N-M+1;--j;s0=s1,s1=*p2++)
    *p0++ = *pM++ ^ (mixBits(s0,s1)>>1) ^ (loBit(s1) ? K : 0U);
  for(pM=state,j=M;--j;s0=s1,s1=*p2++)
    *p0++ = *pM++ ^ (mixBits(s0,s1)>>1) ^ (loBit(s1) ? K : 0U);
  
  s1=state[0],*p0=*pM ^ (mixBits(s0,s1)>>1) ^ (loBit(s1) ? K : 0U);
  s1 ^= (s1 >> 11);
  s1 ^= (s1 << 7) & 0x9D2C5680U;
  s1 ^= (s1 << 15) & 0xEFC60000U;
  return(s1 ^ (s1 >> 18));
}


uint32 randomMT(void)
{
  uint32 y;
  
  if(--left < 0) return(reloadMT());
  
  y = *next++;
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9D2C5680U;
  y ^= (y << 15) & 0xEFC60000U;
  return(y ^ (y >> 18));
}

double ranMT(void)
{
  return (double) randomMT()/NORM;
}

int nearest(double *px,double *py,int k)
{
  int nn=1,i,j;
  double dum=sqrt(2)*L,dum1;
  double dis[9];
  
  for(i=1;i<=S;i++) {
    if(i==k) continue;
    dis[0]=sqrt((px[i]-px[k])*(px[i]-px[k])+(py[i]-py[k])*(py[i]-py[k]));
    dis[1]=sqrt((px[i]-px[k])*(px[i]-px[k])+(py[i]-L-py[k])*(py[i]-L-py[k]));
    dis[2]=sqrt((px[i]-px[k])*(px[i]-px[k])+(py[i]+L-py[k])*(py[i]+L-py[k]));
    dis[3]=sqrt((px[i]-L-px[k])*(px[i]-L-px[k])+(py[i]-py[k])*(py[i]-py[k]));
    dis[4]=sqrt((px[i]-L-px[k])*(px[i]-L-px[k])+(py[i]-L-py[k])*(py[i]-L-py[k]));
    dis[5]=sqrt((px[i]-L-px[k])*(px[i]-L-px[k])+(py[i]+L-py[k])*(py[i]+L-py[k]));
    dis[6]=sqrt((px[i]+L-px[k])*(px[i]+L-px[k])+(py[i]-py[k])*(py[i]-py[k]));
    dis[7]=sqrt((px[i]+L-px[k])*(px[i]+L-px[k])+(py[i]-L-py[k])*(py[i]-L-py[k]));
    dis[8]=sqrt((px[i]+L-px[k])*(px[i]+L-px[k])+(py[i]+L-py[k])*(py[i]+L-py[k]));
    dum1=dis[0]; 
    for(j=1;j<9;j++) {
      if(dis[j] < dum1) dum1 = dis[j];
    }


    if(dum1 < dum) {
      nn=i; dum=dum1;
    }
  }
  
  return(nn);
}

double incremod(double r,double v,double D)
{
  r += v;
  if(r>0) r = fmod(r,D);
  else r = D - fmod(abs(r),D);
  return(r);
}

double orderpar(double *Gam)
{
  int i;
  double p=0;
  
  for(i=1;i<=S;i++) {
    p += cos(2*Pi*Gam[i]);
  }
  p /= S;
  
  return(p);
}

void bhit(double *x,double *y,double *vx,double *vy,double *Gam,double *t)
{
  int i;
  
  /* update phases and positions */
  for(i=1;i<=S;i++) {
	Gam[i] += (*t);
    x[i] += vx[i]*(*t); y[i] += vy[i]*(*t);
	/* Periodical boundary conditions */
	if(x[i] > L) {
		x[i] -= L;
	}
	if(y[i] > L) {
		y[i] -= L;
	}
	if(x[i] < 0) {
		x[i] += L;
	}
	if(y[i] < 0) {
		y[i] += L;
	}    
  }

}

void printOrderpar(FILE *output, double *Gam)
{
	double phase;	
	phase = orderpar(Gam);
	fprintf(output,"%.5g ",phase);
    fprintf(output,"\n");
}