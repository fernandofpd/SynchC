/* swarmsyncRefRan.c - double precision */

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
#include "cokus.h"
#define Pi         2.*asin(1.)

/* ---------------------------------------------------- */

void maxfind(int *p,double *Gam);
void sort(unsigned long n, double a[], int b[]);
/* random number generator */
double ranMT(void);
/* euclid dist p k */
double edist(double *px,double *py,int p,int k);
double cangle(double *px,double *py,double *vx,double *vy,int p,int k);
double orderpar(double *Gam);
int bhit(double *x,double *y,double *vx,double *vy,double *Gam,double *t);


/* Parameter: */
double tau=1.;
/* epsilon, unit number S, box size L, velocity Vel */
double eps=0.1;
int S=20;
double L=400;
double Vel=0.1;
double smin=1.e-6;
double Tmax = 1e7;
double alpha = 18;

int main(int argc,char **argv)
{
	int i,j,k,seed1=7,n=200,p,q;
	int R=1;
	uint32 seed=5;
	double t,phase,phi,dpos,zdum=1,time=0;
	double dposm;
	/* neuron phase Gam, motion angle Phi */ 
	/* xy-position Pxy, xy-velocity Vxy */
	double *Gam,*Phi,*Px,*Py,*Vx,*Vy;
	int *P,*Stime;
	FILE *out1;
	char filename[128]="dat";

	if(opt_flag(&argc,argv,"-h")) {
		printf("\nOptions: -n\n\n");
		exit(0);
	}

	opt_int(&argc,argv,"-n",&n);
	opt_int(&argc,argv,"-N",&S);
	opt_int(&argc,argv,"-s",&seed1);
	opt_double(&argc,argv,"-V",&Vel);
	opt_double(&argc,argv,"-e",&eps);
	opt_double(&argc,argv,"-L",&L);
	opt_double(&argc,argv,"-a",&alpha);
        opt_double(&argc,argv,"-T",&Tmax);
	opt_int(&argc,argv,"-R",&R);
	opt_string(&argc,argv,"-f",filename);	
  
	/* scale max separation and define max angle */
	dposm=2*L/sqrt(S);
	alpha = 3.1416/180*alpha;
  
	out1=fopen(filename,"w"); fclose(out1); out1=fopen(filename,"a");

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
		while(zdum > smin && time < Tmax) {
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
			}
			Gam[p]=0;
      
			for(j=1;j<=S;j++){
				if(j==p) continue;
					dpos=edist(Px,Py,p,j);
					if(dpos <= dposm){
						phi=cangle(Px,Py,Vx,Vy,p,j);
						if(phi <= alpha) {
							Gam[j] *= (1+eps);
							if(Gam[j] >= 1) Gam[j] = 1;
						}
					}
				}
			}
		}
                if(time >= Tmax) time = -1;	// Data is censored
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
  double angle;
  
  /* use scalar product between vel(p) and pos(k)-pos(p) */
  angle = (vx[p]*(px[k]-px[p])+vy[p]*(py[k]-py[p]))/(Vel*edist(px,py,p,k));
  angle = acos(angle); 
 
  return(angle);
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
