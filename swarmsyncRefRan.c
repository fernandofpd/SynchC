/* swarmsyncRefRan.c - double precision */

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <limits.h>

#include "nrutil.h"
#include "vecmatd.h"
#include "randomd.h"
#include "options.h"
#include "cokus.h"
#define Pi         2.*asin(1.)

/* ---------------------------------------------------- */

void maxfind(int *p, double *Gam);
/* random number generator */
double ranMT(void);
void findNeighbors(double *Px, double *Py, double *Vx, double *Vy, int p, int *neigh);
/* find Q nearest neighbors of unit k */
void qNearest(double *Px, double *Py, int p, int *neigh);
/* find neighbors in cone */
void cone(double *Px, double *Py, double *Vx, double *Vy, int p, int *neigh);
/* euclid dist p k */
double edist(double *Px,double *Py,int p,int k);
double cangle(double *Px,double *Py,double *Vx,double *Vy,int p,int k);
/* Calculate Order Parameter */
double orderpar(double *Gam);
int bhit(double *x,double *y,double *Vx,double *Vy,double *Gam,double *t);

/* Parameter: */
double tau = 1.;
/* epsilon, unit number S, box size L, velocity Vel */
double eps = 0.1;
int S = 20;
double L = 400;
double Vel = 0.1;
/* synchronization minimum smin, maximum cycle count Tmax */
double smin = 1.e-6;
double Tmax = 1e7;
char neighborhood[128] = "QNearest";
/* angle alpha and radius r of interaction */
double alpha = 18;
double r = 400;

int Q = 1;

int main(int argc,char **argv)
{
    int i, j, k, seed1 = 7, p, q, R = 1;
    uint32 seed = 5;
    double t, phase, zdum, time;
    /* oscillator phase Gam, motion angle Phi */
    /* xy-position Pxy, xy-velocity Vxy */
    double *Gam, *Phi, *Px, *Py, *Vx, *Vy;
    int *P, *neigh, *Stime;
    FILE *out1;
    char filename[128] = "dat";

    if(opt_flag(&argc,argv,"-h")) {
        printf("\nOptions: -n\n\n");
        exit(0);
    }

    /* Input parameters */
    opt_int(&argc, argv, "-Q", &Q);
    opt_int(&argc, argv, "-N", &S);
    opt_int(&argc, argv, "-R", &R);
    opt_int(&argc, argv, "-s", &seed1);
    opt_double(&argc, argv, "-V", &Vel);
    opt_double(&argc, argv, "-e", &eps);
    opt_double(&argc, argv, "-L", &L);
    opt_double(&argc,argv,"-a",&alpha);
    opt_double(&argc, argv, "-r", &r);
    opt_double(&argc,argv,"-T",&Tmax);
    opt_string(&argc, argv, "-f", filename);
    opt_string(&argc, argv, "-n", neighborhood);

    alpha = Pi/180*alpha;

    out1 = fopen(filename, "w"); fclose(out1); out1 = fopen(filename, "a");
 
    Gam = dvector(1, S);
    Phi = dvector(1, S);
    Px = dvector(1, S);
    Py = dvector(1, S);
    Vx = dvector(1, S);
    Vy = dvector(1, S);
    Stime = ivector(1, R);
    neigh = ivector(1, S);

    /* initialize random numbers */
    seed = (uint32)(seed1); seedMT(seed);

    /* Loop over the Runs */
    for(j = 1; j <= R; j++) {
        zdum = 1; time = 0;

        /* random initial values */
        for(i = 1; i <= S; i++) {
            Gam[i] = ranMT(); Phi[i] = 2*Pi*ranMT();
            Px[i] = L*ranMT(); Py[i] = L*ranMT();
            Vx[i] = Vel*cos(Phi[i]); Vy[i] = Vel*sin(Phi[i]);
        }

        while(zdum > smin && time < Tmax) {
            q = 0;
            /* index of next firing unit p, time t */
            maxfind(&p, Gam); t = 1 - Gam[p];
            /* update position and velocities */
            /* either next unit fires (q=1) or hits the wall (q=0) */
            q = bhit(Px, Py, Vx, Vy, Gam, &t);

            if(q > 0) { 
                /* when reference unit 1 fires find out order parameter */
                if(p == 1) {
                    phase = orderpar(Gam);
                    zdum = 1 - phase;
                    time++;
                }
                Gam[p] = 0;            
                
                findNeighbors(Px, Py, Vx, Vy, p, neigh);            
                for(k = 1; k <= S; k++) {
                    if(neigh[k] == 1) {
                        Gam[k] *= (1 + eps);
                        if(Gam[k] >= 1) Gam[k] = 1;
                    }
                }
            }
        }
        if(time >= Tmax) time = -1;    // Data is censored
        Stime[j] = time;
    }
    /* Save results to file */
    for(i = 1; i <= R; i++) fprintf(out1, "%d\n", Stime[i]);
    
    /* Free memory */
    fclose(out1);
    free_dvector(Gam, 1, S);
    free_dvector(Phi, 1, S);
    free_dvector(Px, 1, S);
    free_dvector(Py, 1, S);
    free_dvector(Vx, 1, S);
    free_dvector(Vy, 1, S);
    free_ivector(Stime, 1, R);
    free_ivector(neigh, 1, S);
    free_ivector(P,1,S);

return 0;
}

void maxfind(int *p, double *Gam)
{
    int i;
    *p=1;
    for(i = 2; i <= S; i++) if(Gam[i] > Gam[*p]) *p = i;
}

double ranMT(void)
{
    return (double) randomMT()/NORM;
}

void findNeighbors(double *Px, double *Py, double *Vx, double *Vy, int p, int *neigh) {
    int i;
    // Initialize neighborhood to zero
    for(i = 1; i <= S; i++) neigh[i] = 0;
    if(strcmp(neighborhood,"QNearest") == 0) qNearest(Px, Py, p, neigh);
    else if(strcmp(neighborhood,"ConeOut") == 0) cone(Px, Py, Vx, Vy, p, neigh);
}

void qNearest(double *Px, double *Py, int p, int *neigh) {
    int i, j, foo;
    double dum = 2*L*L, dum1;
    double dis[9];
    double *dis2all;
    dis2all = dvector(1, S);
    dis2all[p] = dum;

    /* Calculate distances to all other oscillators */
    for(i = 1; i <= S; i++) {
        if(i == p) continue;
        /* Take into account the periodic boundary conditions */
        dis[0] = (Px[i]-Px[p])*(Px[i]-Px[p])+(Py[i]-Py[p])*(Py[i]-Py[p]);
        dis[1] = (Px[i]-Px[p])*(Px[i]-Px[p])+(Py[i]-L-Py[p])*(Py[i]-L-Py[p]);
        dis[2] = (Px[i]-Px[p])*(Px[i]-Px[p])+(Py[i]+L-Py[p])*(Py[i]+L-Py[p]);
        dis[3] = (Px[i]-L-Px[p])*(Px[i]-L-Px[p])+(Py[i]-Py[p])*(Py[i]-Py[p]);
        dis[4] = (Px[i]-L-Px[p])*(Px[i]-L-Px[p])+(Py[i]-L-Py[p])*(Py[i]-L-Py[p]);
        dis[5] = (Px[i]-L-Px[p])*(Px[i]-L-Px[p])+(Py[i]+L-Py[p])*(Py[i]+L-Py[p]);
        dis[6] = (Px[i]+L-Px[p])*(Px[i]+L-Px[p])+(Py[i]-Py[p])*(Py[i]-Py[p]);
        dis[7] = (Px[i]+L-Px[p])*(Px[i]+L-Px[p])+(Py[i]-L-Py[p])*(Py[i]-L-Py[p]);
        dis[8] = (Px[i]+L-Px[p])*(Px[i]+L-Px[p])+(Py[i]+L-Py[p])*(Py[i]+L-Py[p]);
        dum1 = dis[0]; 
        for(j = 1; j < 9; j++) {
            if(dis[j] < dum1) dum1 = dis[j];
        }
        dis2all[i] = dum1;
    }

    /* Find Q nearest neighbors */
    for(i = 1; i <= Q; i++) {
        dum1 = dum;
        for(j = 1; j <= S; j++) {
            if(dis2all[j] < dum1) {
                foo = j;
                dum1 = dis2all[j];
            }
        }
        dis2all[foo] = dum;
        neigh[foo] = 1;
    }
    free_dvector(dis2all, 1, S);
}

void cone(double *Px, double *Py, double *Vx, double *Vy, int p, int *neigh) {
    int j;
    double phi, dpos;
    for(j = 1; j <= S; j++) {
        if(j == p) continue;
        dpos = edist(Px, Py, p, j);
        if(dpos <= r) {
            phi = cangle(Px, Py, Vx, Vy, p, j);
            if(phi <= alpha) neigh[j] = 1;            
        }
    }
}

double edist(double *Px,double *Py,int p,int k)
{
  double dis;
  dis = sqrt((Px[p]-Px[k])*(Px[p]-Px[k])+(Py[p]-Py[k])*(Py[p]-Py[k]));
  return(dis);
}

double cangle(double *Px,double *Py,double *Vx,double *Vy,int p,int k)
{
  double angle;
  
  /* use scalar product between vel(p) and pos(k)-pos(p) */
  angle = (Vx[p]*(Px[k]-Px[p])+Vy[p]*(Py[k]-Py[p]))/(Vel*edist(Px,Py,p,k));
  angle = acos(angle); 
 
  return(angle);
}

double orderpar(double *Gam)
{
    int i;
    double p = 0;
    for(i = 1; i <= S; i++) p += cos(2*Pi*Gam[i]);    
    p /= S;
    return(p);
}

int bhit(double *x,double *y,double *Vx,double *Vy,double *Gam,double *t)
{
    int i, q = 1, qi = 1;
    double tb = 2, dum;
    
    /* compute time for next boundary hit */
    for(i = 1; i <= S; i++) {
        if(Vx[i] > 0) dum = (L-x[i])/Vx[i];
        if(Vx[i] < 0) dum = (-x[i])/Vx[i];
        if(Vx[i] == 0) dum = 2;
        if(dum < tb) {
            tb = dum; q = i;
        }
    }
    for(i = S+1; i <= 2*S; i++) {
        if(Vy[i-S] > 0) dum = (L-y[i-S])/Vy[i-S];
        if(Vy[i-S] < 0) dum = (-y[i-S])/Vy[i-S];
        if(Vy[i-S] == 0) dum = 2;
        if(dum < tb) {
            tb = dum; q = i;
        }
    }
    /* boundary hit before next firing? */
    if(tb <= *t) {
        qi = 0; *t = tb;
    }
    /* update phases and positions */
    for(i=1;i<=S;i++) {
        x[i] += Vx[i]*(*t); y[i] += Vy[i]*(*t);
        Gam[i] += (*t);
    }
    /* if wall is hit do random reflection */
    /* distinguish x-wall (q <= S) and y-wall */
    if(qi == 0) {
        dum = Pi*ranMT();
        if(q <= S) {
            Vx[q] = -(Vx[q]>=0)*sin(dum)*Vel + (Vx[q]<0)*sin(dum)*Vel;
            Vy[q] = Vel*cos(2*dum);
        } 
        else {
            Vy[q-S] = -(Vy[q-S]>=0)*sin(dum)*Vel + (Vy[q-S]<0)*sin(dum)*Vel;
            Vx[q-S] = Vel*cos(2*dum);
        }
    }
    return(qi);
}
