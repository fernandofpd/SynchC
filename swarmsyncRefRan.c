/* swarmsyncRefRan.c - double precision */
/* unit free */
/* random re-orientation at boundary */
/* interaction with neurons in cone */
/* Verwendet Routinen aus Numerical Recipes */
/* Ruediger Zillmer, May 2013 */

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
/* euclid dist p k */
double edist(double *px,double *py,int p,int k);
double cangle(double *px,double *py,double *vx,double *vy,int p,int k);
double orderpar(double *Gam);
int bhit(double *x,double *y,double *vx,double *vy,double *Gam,double *t);


/* Parameter: */
double tau=1.;
double dphi=0.1;
/* epsilon, unit number S, box size L, velocity Vel */
double eps=0.1;
int S=20;
double L=400;
double Vel=0.1;
double smin=1.e-6;

int main(int argc,char **argv)
{
  int i,j,k,seed1=7,n=200,p,q;
  int R=1;
  uint32 seed=5;
  double t,phase,cphi,dpos,zdum=1,time=0;
  double dposm,cphim,alpha=18;
  /* neuron phase Gam, motion angle Phi */ 
  /* xy-position Pxy, xy-velocity Vxy */
  double *Gam,*Phi,*Px,*Py,*Vx,*Vy;
  int *P,*Stime;
  FILE *out1;

  if(opt_flag(&argc,argv,"-h")) {
    printf("\nOptions: -n\n\n");
    exit(0);
  }

  /*sprintf(file1,"map");
  opt_string(&argc,argv,"-o",file1);
  strcat(file1,".tim");

  opt_save(&argc,argv,file1,"w"); */

  opt_int(&argc,argv,"-n",&n);
  opt_int(&argc,argv,"-N",&S);
  opt_int(&argc,argv,"-s",&seed1);
  opt_double(&argc,argv,"-dp",&dphi);
  opt_double(&argc,argv,"-V",&Vel);
  opt_double(&argc,argv,"-e",&eps);
  opt_double(&argc,argv,"-L",&L);
  opt_double(&argc,argv,"-a",&alpha);
  opt_int(&argc,argv,"-R",&R);
  
  /* scale max separation and define max angle */
  dposm=2*L/sqrt(S);
  cphim=cos(3.1416/180*alpha);
  
  out1=fopen("dat","w"); fclose(out1); out1=fopen("dat","a");

  /*TT = ivector(1,S);*/
  P = ivector(1,S);
  Gam = dvector(1,S);
  Phi = dvector(1,S);
  Px = dvector(1,S);
  Py = dvector(1,S);
  Vx = dvector(1,S);
  Vy = dvector(1,S);
  Stime = ivector(1,R);
  
  /* initialize random numbers */
  seed=(uint32)(seed1); seedMT(seed);
  
  for(k=1;k<=R;k++) {
    zdum=1; time=0;
  
  /* random initial values */
  for(i=1;i<=S;i++) {
    Gam[i]=ranMT(); Phi[i]=2*Pi*ranMT();
    Px[i]=L*ranMT(); Py[i]=L*ranMT();
    Vx[i]=Vel*cos(Phi[i]); Vy[i]=Vel*sin(Phi[i]);
  } 
  
  /*for(i=1;i<=n;i++) {*/
  while(zdum > smin) {
    q=0;
    /* index of next firing unit p, time t */
    maxfind(&p,Gam); t = 1-Gam[p];
    /* update position and velocities */
    /* either next unit fires (q=1) or hits the wall (q=0) */
    q = bhit(Px,Py,Vx,Vy,Gam,&t);
            
    if(q>0) { 
      /* when reference unit 1 fires print out order parameter */
      if(p==1) {
        phase = orderpar(Gam);
	zdum = 1-phase; time++;
        //for(j=1;j<=S;j++) fprintf(out1,"%.5g %.5g %.5g ",Px[j],Vx[j],Gam[j]);
        //fprintf(out1,"%.5g %.5g ",Px[1],Py[1]);
        //fprintf(out1,"%.5g ",phase);
        //fprintf(out1,"\n");
      }
      Gam[p]=0;
      
      for(j=1;j<=S;j++){
        if(j==p) continue;
        dpos=edist(Px,Py,p,j);
	if(dpos <= dposm){
	  cphi=cangle(Px,Py,Vx,Vy,p,j);
	  if(cphi >= cphim) {
	    Gam[j] *= (1+eps);
            if(Gam[j] >= 1) Gam[j] = 1;
          }
      	}
      }
    }
  }
  Stime[k] = time;
  }
  
  for(i=1;i<=R;i++) fprintf(out1,"%d\n",Stime[i]); 
  
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

double edist(double *px,double *py,int p,int k)
{
  double dis;
  
  dis=sqrt((px[p]-px[k])*(px[p]-px[k])+(py[p]-py[k])*(py[p]-py[k]));
  
  return(dis);
}

double cangle(double *px,double *py,double *vx,double *vy,int p,int k)
{
  double dis;
  
  dis = sqrt((px[p]-px[k])*(px[p]-px[k])+(py[p]-py[k])*(py[p]-py[k]));
  /* use scalar product between vel(p) and pos(k)-pos(p) */
  dis = (vx[p]*(px[k]-px[p])+vy[p]*(py[k]-py[p]))/(Vel*dis);
  
  return(dis);
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

int bhit(double *x,double *y,double *vx,double *vy,double *Gam,double *t)
{
  int i,q=1,qi=1;
  double tb=2,dum;
  
  /* compute time for next boundary hit */
  for(i=1;i<=S;i++) {
    if(vx[i] > 0) dum = (L-x[i])/vx[i];
    if(vx[i] < 0) dum = (-x[i])/vx[i];
    if(vx[i] == 0) dum = 2;
    if(dum < tb) {
      tb = dum; q = i;
    }
  }
  for(i=(S+1);i<=2*S;i++) {
    if(vy[i-S] > 0) dum = (L-y[i-S])/vy[i-S];
    if(vy[i-S] < 0) dum = (-y[i-S])/vy[i-S];
    if(vy[i-S] == 0) dum = 2;
    if(dum < tb) {
      tb = dum; q = i;
    }
  }
  /* boundary hit before next firing? */
  if(tb <= *t) {
    qi=0; *t=tb;
  }
  /* update phases and positions */
  for(i=1;i<=S;i++) {
    x[i] += vx[i]*(*t); y[i] += vy[i]*(*t);
    Gam[i] += (*t);
  }
  /* if wall is hit do random reflection */
  /* distinguish x-wall (q <= S) and y-wall */
  if(qi == 0) {
    dum = Pi*ranMT();
    if(q <= S) {
      vx[q] = -(vx[q]>=0)*sin(dum)*Vel + (vx[q]<0)*sin(dum)*Vel;
      vy[q] = Vel*cos(2*dum);
    } 
    else {
      vy[q-S] = -(vy[q-S]>=0)*sin(dum)*Vel + (vy[q-S]<0)*sin(dum)*Vel;
      vx[q-S] = Vel*cos(2*dum);
    }
  }
    
  return(qi);
}