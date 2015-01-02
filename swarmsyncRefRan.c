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
void findNeighbors(double *px, double *py, int k, int Q, int *neigh);
/* find Q nearest neighbors of unit k */
void Qnearest(double *px, double *py, int k, int Q, int *neigh);
double orderpar(double *Gam);
void bhit(double *x, double *y, double *vx, double *vy, double *Gam, double *t);

/* Parameter: */
double tau = 1.;
double dphi = 0.1;
/* epsilon, unit number S, box size L, velocity Vel */
double eps = 0.1;
int S = 20;
double L = 400;
double Vel = 0.1;
double smin = 1.e-6;
int maxT = 100000;
char neighborhood[128] = "QNearest";

int main(int argc,char **argv)
{
    int i, j, k, seed1 = 7, p, Q = 1, R = 1, time = 0;
    uint32 seed = 5;
    double t, phase, zdum;
    /* neuron phase Gam, motion angle Phi */
    /* xy-position Pxy, xy-velocity Vxy */
    double *Gam, *Phi, *Px, *Py, *Vx, *Vy;
    int *neigh, *Stime;
    FILE *out1;
    char filename[128] = "dat";
    if(INT_MAX < maxT) maxT = INT_MAX -1;

    /* Input parameters */
    opt_int(&argc, argv, "-Q", &Q);
    opt_int(&argc, argv, "-N", &S);
    opt_int(&argc, argv, "-R", &R);
    opt_int(&argc, argv, "-s", &seed1);
    opt_double(&argc, argv, "-dp", &dphi);
    opt_double(&argc, argv, "-V", &Vel);
    opt_double(&argc, argv, "-e", &eps);
    opt_double(&argc, argv, "-L", &L);    
    opt_string(&argc, argv, "-f", filename);
    opt_string(&argc, argv, "-n", neighborhood);

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
    seed=(uint32)(seed1); seedMT(seed);
	
    for(j=1;j<=R;j++) {
        zdum=1; time=0;

        /* random initial values */
        for(i = 1; i <= S; i++) {
            Gam[i] = ranMT(); Phi[i] = 2*Pi*ranMT();
            Px[i] = L*ranMT(); Py[i] = L*ranMT();
            Vx[i] = Vel*cos(Phi[i]); Vy[i] = Vel*sin(Phi[i]);
        } 

        while(zdum > smin && time < maxT) {
            /* index of next firing unit p, time t */
            maxfind(&p, Gam); t = 1-Gam[p];
            /* update position and velocities */
            bhit(Px, Py, Vx, Vy, Gam, &t);

            /* when reference unit 1 fires find out order parameter */
            if(p == 1) {
                phase = orderpar(Gam);
                zdum = 1-phase;
                time++;
            }
            Gam[p] = 0;            

            findNeighbors(Px, Py, p, Q, neigh);            
            for(k = 1; k <= S; k++) {
                if(neigh[k] == 1) {
                    Gam[k] *= (1+eps);
                    if(Gam[k] >= 1) Gam[k] = 1;
                }
            }
        }
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

void findNeighbors(double *px, double *py, int k, int Q, int *neigh) {
    int i;
    // Initialize neighborhood to zero
    for(i = 1; i <= S; i++) neigh[i] = 0;
    if(strcmp(neighborhood,"QNearest") == 0) Qnearest(px, py, k, Q, neigh);
}

void Qnearest(double *px, double *py, int k, int Q, int *neigh) {
    int i, j, foo;
    double dum = sqrt(2.)*L, dum1;
    double dis[9];
    double *dis2all;
    dis2all = dvector(1, S);
    dis2all[k] = dum;

    /* Calculate distances to all other oscillators */
    for(i=1;i<=S;i++) {
        if(i==k) continue;
        /* Take into account the periodic boundary conditions */
        dis[0] = sqrt((px[i]-px[k])*(px[i]-px[k])+(py[i]-py[k])*(py[i]-py[k]));
        dis[1] = sqrt((px[i]-px[k])*(px[i]-px[k])+(py[i]-L-py[k])*(py[i]-L-py[k]));
        dis[2] = sqrt((px[i]-px[k])*(px[i]-px[k])+(py[i]+L-py[k])*(py[i]+L-py[k]));
        dis[3] = sqrt((px[i]-L-px[k])*(px[i]-L-px[k])+(py[i]-py[k])*(py[i]-py[k]));
        dis[4] = sqrt((px[i]-L-px[k])*(px[i]-L-px[k])+(py[i]-L-py[k])*(py[i]-L-py[k]));
        dis[5] = sqrt((px[i]-L-px[k])*(px[i]-L-px[k])+(py[i]+L-py[k])*(py[i]+L-py[k]));
        dis[6] = sqrt((px[i]+L-px[k])*(px[i]+L-px[k])+(py[i]-py[k])*(py[i]-py[k]));
        dis[7] = sqrt((px[i]+L-px[k])*(px[i]+L-px[k])+(py[i]-L-py[k])*(py[i]-L-py[k]));
        dis[8] = sqrt((px[i]+L-px[k])*(px[i]+L-px[k])+(py[i]+L-py[k])*(py[i]+L-py[k]));
        dum1=dis[0]; 
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

double orderpar(double *Gam)
{
    int i;
    double p = 0;
    for(i = 1; i <= S; i++) p += cos(2*Pi*Gam[i]);    
    p /= S;
    return(p);
}

void bhit(double *x,double *y,double *vx,double *vy,double *Gam,double *t)
{
    int i;  
    /* update phases and positions */
    for(i = 1; i <= S; i++) {
        Gam[i] += (*t);
        x[i] += vx[i]*(*t); y[i] += vy[i]*(*t);
        /* Periodical boundary conditions */
        if(x[i] > L) x[i] -= L;        
        if(y[i] > L) y[i] -= L;
        if(x[i] < 0) x[i] += L;
        if(y[i] < 0) y[i] += L;    
    }
}
