/* ---------- swarmSync.c --------- */
/* Swarm Synchronization of S particles moving with speed vel in a square environment of length L */
/* Particles are linear integrate and fire oscillators with period tau */
/* When an oscillator's  phase reaches 1 it is reset to 0 and interacts with its neighbors */
/* This interaction is a multiplicative update of the neighbor's phase by a factor (1 + epsilon) */
/* Several type of neighborhoods are included */
/* The system is updated either upon a wall-hit or a firing event. It is a pure continuous time simulation */
/* Fernando Perez-Diaz and Ruediger Zillmer */

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <limits.h>

#include "nrutil.h"
#include "options.h"
#include "cokus.h"
#define Pi 2.*asin(1.)

/* ----------- Function Declarations ----------- */

void maxfind(int *p, double *phase);
/* Random number generator */
double ranMT(void);
/* Find Neighbors */
void findNeighbors(double *Px, double *Py, double *Vx, double *Vy, int p, int *neigh);
/* Find Q nearest neighbors of unit k */
void qNearest(double *Px, double *Py, int p, int *neigh);
/* Find neighbors in cone */
void cone(double *Px, double *Py, double *Vx, double *Vy, int p, int *neigh, char *direction);
/* Distance and angle between p & k */
double edist(double *Px, double *Py, int p, int k);
double cangle(double *Px, double *Py, double *Vx, double *Vy, int p, int k);
/* Calculate Order Parameter */
double orderpar(double *phase);
/* Calculate next firing or wall hit */
int bhit(double *x, double *y, double *Vx, double *Vy, double *phase, double *t);

/* ----------  Parameters ----------- */

/* Tau, epsilon, unit number S, box size L, velocity Vel */
double tau = 1.;
double eps = 0.1;
int S = 20;
double L = 400;
double Vel = 0.1;
/* Synchronization minimum smin, maximum cycle count Tmax */
double smin = 1.e-6;
double Tmax = 1e7;
/* Type of neighborhood */
char neighborhood[128] = "QNearest";
/* Angle alpha and radius r of interaction */
double alpha = 18;
double r = 400;
/* Number of nearest neighbors */
int Q = 1;

/* ---------- Main ----------- */

int main(int argc,char **argv)
{
    /* ------ Variable Declarations -----*/
    int i, j, k, seed1 = 7, p, q, R = 1, calcConn;
    uint32 seed = 5;
    double t, ordPar, zdum, time, firingTime;
    /* Oscillator phase, motion angle phi */
    /* xy-position Pxy, xy-velocity Vxy */
    double *phase, *phi, *Px, *Py, *Vx, *Vy;
    int *neigh, *Stime, **conn;
    FILE *out1, *out2;
    char filename[128] = "dat", filename2[128] = "conn";

    /* ----- Input parameters ----- */
    
    calcConn = opt_flag(&argc, argv,"-c"); 
    opt_int(&argc, argv, "-Q", &Q);
    opt_int(&argc, argv, "-N", &S);
    opt_int(&argc, argv, "-R", &R);
    opt_int(&argc, argv, "-s", &seed1);
    opt_double(&argc, argv, "-V", &Vel);
    opt_double(&argc, argv, "-e", &eps);
    opt_double(&argc, argv, "-L", &L);
    opt_double(&argc, argv,"-a", &alpha);
    opt_double(&argc, argv, "-r", &r);
    opt_double(&argc,argv,"-T",&Tmax);
    opt_string(&argc, argv, "-f", filename);
    opt_string(&argc, argv, "-n", neighborhood);

    /* ----- Initializations ----- */

    alpha = Pi/180*alpha;

    out1 = fopen(filename, "w"); fclose(out1); out1 = fopen(filename, "a");
 
    phase = dvector(1, S);
    phi = dvector(1, S);
    Px = dvector(1, S);
    Py = dvector(1, S);
    Vx = dvector(1, S);
    Vy = dvector(1, S);
    Stime = ivector(1, R);
    neigh = ivector(1, S);
    conn = imatrix(1, S, 1, S);

    /* Initialize random numbers */
    seed = (uint32)(seed1); seedMT(seed);

    /* Loop over the Runs */
    for(j = 1; j <= R; j++) {
        if(calcConn) {
            sprintf(filename2, "%sconn%d", filename, j);
            out2 = fopen(filename2, "w"); fclose(out2); out2 = fopen(filename2, "a");
        }

        zdum = 1; time = 0;
        firingTime = 0;

        /* Random initial values */
        for(i = 1; i <= S; i++) {
            phase[i] = ranMT(); phi[i] = 2*Pi*ranMT();
            Px[i] = L*ranMT(); Py[i] = L*ranMT();
            Vx[i] = Vel*cos(phi[i]); Vy[i] = Vel*sin(phi[i]);
        }

        /* While not yet synchronized or before Tmax cycles have elapsed */
        while(zdum > smin && time < Tmax) {
            q = 0;
            /* Index of next firing unit p, time t */
            maxfind(&p, phase); t = 1 - phase[p];
            /* Update position and velocities */
            /* either next unit fires (q=1) or hits the wall (q=0) */
            q = bhit(Px, Py, Vx, Vy, phase, &t);

            /* If unit fires */
            if(q > 0) {
                if(t > 0 && calcConn == 1) {
                    /* Save connectivity matrix */
                    if(firingTime > 0) {
                        fprintf(out2, "%f\t", firingTime);
                        for(i = 1; i <= S; i++) {
                            for(k = 1; k <= S; k++) {
                               fprintf(out2, "%d\t", conn[i][k]);
                            }
                        }
                        fprintf(out2, "\n");
                    }
                    /* Reset matrix to zero */                    
                    for(i = 1; i <= S; i++) {
                        for(k = 1; k <= S; k++) {
                            conn[i][k] = 0;
                        }
                    }
                    firingTime += t; 
                }
                /* When reference unit 1 fires calculate order parameter */
                if(p == 1) {
                    ordPar = orderpar(phase);
                    zdum = 1 - ordPar;
                    time++;
                }
                /* Reset phase */
                phase[p] = 0; 
                /* Find neighbors and update their phase */
                findNeighbors(Px, Py, Vx, Vy, p, neigh); 
                for(k = 1; k <= S; k++) {
                    if(neigh[k] == 1) {
                        phase[k] *= (1 + eps);
                        if(phase[k] >= 1) phase[k] = 1;
                        /* Update connectivity matrix */
                        if(calcConn == 1) conn[p][k] += 1;
                    }
                }
            }
        }
        /* If not synchronized in Tmax cycles, censor data */
        if(time >= Tmax) time = -1;
        /* Save the Synchronization times */
        Stime[j] = time;
        fprintf(out1, "%d\n", Stime[j]);

        if(calcConn) fclose(out2);
    }
 
    /* ----- Free memory ----- */
    fclose(out1);
    free_dvector(phase, 1, S);
    free_dvector(phi, 1, S);
    free_dvector(Px, 1, S);
    free_dvector(Py, 1, S);
    free_dvector(Vx, 1, S);
    free_dvector(Vy, 1, S);
    free_ivector(Stime, 1, R);
    free_ivector(neigh, 1, S);
    free_imatrix(conn, 1, S, 1, S);

    return 0;
}

/* ----- Find index p corresponding to maximum phase ----- */
void maxfind(int *p, double *phase)
{
    int i;
    *p = 1;
    for(i = 2; i <= S; i++) if(phase[i] > phase[*p]) *p = i;
}

/* ----- Random number generator ----- */
double ranMT(void)
{
    return (double) randomMT()/NORM;
}

/* ------ Find all neighbors ----- */
void findNeighbors(double *Px, double *Py, double *Vx, double *Vy, int p, int *neigh) {
    int i;
    /* Initialize neighborhood to zero */
    for(i = 1; i <= S; i++) neigh[i] = 0;
    /* Use appropriate function according to neighborhood */
    if(strcmp(neighborhood, "QNearest") == 0) qNearest(Px, Py, p, neigh);
    else if(strcmp(neighborhood, "ConeOut") == 0) cone(Px, Py, Vx, Vy, p, neigh, "Out");
    else if(strcmp(neighborhood, "ConeIn") == 0) cone(Px, Py, Vx, Vy, p, neigh, "In");
}

/* ----- Find Q nearest neighbors ----- */
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
        dis2all[i] = edist(Px, Py, p, i);
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

/* ----- Find Neighbors in a Cone ---- */
/* If direction == "Out" : p's neighbors will be the agents that lie inside its cone -- Cone of Influence */
/* If direction == "In" : p's neighbors will be the agents that see p inside their cone -- Cone of Vision */
void cone(double *Px, double *Py, double *Vx, double *Vy, int p, int *neigh, char *direction) {
    int j;
    double phi, dpos;
    for(j = 1; j <= S; j++) {
        if(j == p) continue;
        dpos = edist(Px, Py, p, j);
        if(dpos <= r) {
            if (direction == "Out") phi = cangle(Px, Py, Vx, Vy, p, j);
            if (direction == "In") phi = cangle(Px, Py, Vx, Vy, j, p);
            if(phi <= alpha) neigh[j] = 1;            
        }
    }
}

/* ----- Euclidian distance between p & k ----- */
double edist(double *Px, double *Py, int p, int k)
{
    double dis;
    dis = sqrt((Px[p]-Px[k])*(Px[p]-Px[k]) + (Py[p]-Py[k])*(Py[p]-Py[k]));
    return(dis);
}

/* ----- Angle between p->k and p's direction of motion ----- */
double cangle(double *Px,double *Py,double *Vx,double *Vy,int p,int k)
{
    double angle;

    /* Use scalar product between vel(p) and pos(k)-pos(p) */
    angle = (Vx[p]*(Px[k]-Px[p])+Vy[p]*(Py[k]-Py[p]))/(Vel*edist(Px,Py,p,k));
    angle = acos(angle); 
 
    return(angle);
}

/* ----- Calculate Order Parameter ----- */
double orderpar(double *phase)
{
    int i;
    double p = 0;
    for(i = 1; i <= S; i++) p += cos(2*Pi*phase[i]);    
    p /= S;
    return(p);
}

/* ----- Calculate next firing or wall hit ---- */
/* Update position and velocities */
/* either next unit fires (q=1) or hits the wall (q=0) */
int bhit(double *x, double *y, double *Vx, double *Vy, double *phase, double *t)
{
    int i, q = 1, qi = 1;
    double tb = 2, dum;
    
    /* Compute time for next boundary hit */
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
    /* Boundary hit before next firing? */
    if(tb <= *t) {
        qi = 0; *t = tb;
    }
    /* Update phases and positions */
    for(i = 1; i <= S; i++) {
        x[i] += Vx[i]*(*t); y[i] += Vy[i]*(*t);
        phase[i] += (*t);
    }
    /* If wall is hit do random reflection */
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
