/* ---------- swarmSync.c --------- */
/* Swarm Synchronization of S particles moving with speed "speed" in a square environment of length "L" */
/* Particles are linear integrate and fire oscillators with period "tau" */
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
#define PI 2.*asin(1.)

/* ----------- Function Declarations ----------- */

void maxfind(int *p, double *phase);
/* Random number generator */
double ranMT(void);
/* ----- Random Normal number generator ----- */
double ranNorm(void);
/* N-Dimensional directional cosines */
void cosines(double **theta, int i);
/* Find Neighbors */
void findNeighbors(double **pos, double **vel, int p, int *neigh);
/* Find Q nearest neighbors of unit k */
void qNearest(double **pos, int p, int *neigh);
/* Find neighbors in cone */
void cone(double **pos, double **vel, int p, int *neigh, char *direction);
/* Distance and angle between p & k */
double edist(double **pos, int p, int k);
double cangle(double **pos, double **vel,int p,int k);
/* Calculate Order Parameter */
double orderpar(double *phase);
/* Update positions until next firing */
void unboundedMove(double **pos, double **vel, double t);
void boundedMove(double **pos, double **vel, double **theta, double t);

/* ----------  Parameters ----------- */

/* Tau, epsilon, unit number S, box size L, velocity speed */
double tau = 1.;
double eps = 0.1;
int S = 20;
double L = 400;
double speed = 0.1;
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

int Dims = 2;
int boundary;

/* ---------- Main ----------- */

int main(int argc,char **argv)
{
    /* ------ Variable Declarations -----*/
    int i, j, k, d, seed1 = 7, p, R = 1, calcConn, reorientAtInteraction, reorientAtFiring;
    uint32 seed = 5;
    double t, ordPar, zdum, time, firingTime;
    /* Oscillator phase, motion angle theta */
    /* xy-position Pxy, xy-velocity Vxy */
    double *phase, **theta, **pos,  **vel;
    int *neigh, *Stime;
    FILE *out1, *out2;
    char filename[128] = "dat", filename2[128] = "conn";

    /* ------ Variable Declarations -----*/
    /* ----- Input parameters ----- */
    
    calcConn = opt_flag(&argc, argv, 2, "-c", "--conn");                               // Calculate connectivity (if flag is on)
    boundary = opt_flag(&argc, argv, 2, "-b", "--bounded");                            // Bounded environment (bounce at walls) (if flag is on) or has periodical boundary conditions
    reorientAtInteraction = opt_flag(&argc, argv, 1, "--reorientAtInteraction");       // Reorient upon receiving an interaction 
    reorientAtFiring = opt_flag(&argc, argv, 1, "--reorientAtFiring");                 // Reorient upon firing
    opt_int(&argc, argv, &Q, 2, "-Q", "--qnearest");                                   // Number of Nearest Neighbors
    opt_int(&argc, argv, &S, 2, "-N", "--agents");                                     // Total Number of Agents
    opt_int(&argc, argv, &R, 2, "-R", "--runs");                                       // Number of Runs    
    opt_int(&argc, argv, &Dims, 2, "-D", "--dimensions");                              // Number of Runs    
    opt_int(&argc, argv, &seed1, 2, "-s", "--seed");                                   // Random Number Seed    
    opt_double(&argc, argv, &speed, 2, "-V", "--speed");                               // Speed of agents
    opt_double(&argc, argv, &eps, 2, "-e", "--epsilon");                               // Interaction parameter Epsilon
    opt_double(&argc, argv, &L, 2, "-L", "--length");                                  // Length of the environment
    opt_double(&argc, argv, &alpha, 2, "-a", "--alpha");                               // Angle of Interaction (In/Out)
    opt_double(&argc, argv, &r, 2, "-r", "--radius");                                  // Radius of Interaction (In/Out)
    opt_double(&argc, argv, &Tmax, 2, "-T", "--tmax");                                 // Censoring threshold (if the system is not synchronized before Tmax stop runing)
    opt_string(&argc, argv, filename, 2, "-f", "--filename");                          // Output filename
    opt_string(&argc, argv, neighborhood, 2, "-n", "--neighborhood");                  // Define the neighborhood model: {QNearest,ConeOut,ConeIn}

    /* ----- Initializations ----- */

    alpha = PI/180*alpha;

    out1 = fopen(filename, "w"); fclose(out1); out1 = fopen(filename, "a");
 
    phase = dvector(1, S);
    theta = dmatrix(1, S, 1, Dims);
    pos = dmatrix(1, S, 1, Dims);
    vel = dmatrix(1, S, 1, Dims);
    Stime = ivector(1, R);
    neigh = ivector(1, S);

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
            phase[i] = ranMT(); 
            cosines(theta, i);
            for(d = 1; d <= Dims; d++) {
                pos[i][d] = L*ranMT();
                vel[i][d] = speed*theta[i][d];
            } 
        }

        /***** While not yet synchronized or before Tmax cycles have elapsed *****/
        while(zdum > smin && time < Tmax) {
            /***** Update phases and perform movement *****/
            maxfind(&p, phase); t = 1 - phase[p];       // Index of next firing unit p, time t
            for(i = 1; i <= S; i++) phase[i] += t;      // Update phases
            // Update positions until next firing
            if(t > 0) {
                if(boundary) boundedMove(pos, vel, theta, t);
                else unboundedMove(pos, vel, t);
            }

            /***** Firings and Interactions with Neighbors *****/
            /* When reference unit 1 fires calculate order parameter */
            if(p == 1) {
                ordPar = orderpar(phase);
                zdum = 1 - ordPar;
                time++;
            }            
            phase[p] = 0;   // Reset phase
            /* If flag is active reorient upon firing */
            if(reorientAtFiring) {
                cosines(theta, p);
                for(d = 1; d <= Dims; d++) vel[p][d] = speed*theta[p][d]; 
            }
            /* Find neighbors and update their phase */
            findNeighbors(pos, vel, p, neigh); 
            for(k = 1; k <= S; k++) {
                if(neigh[k] == 1) {  
                    phase[k] *= (1 + eps);
                    if(phase[k] >= 1) phase[k] = 1;
                    /* If flag is active reorient upon receiving an interaction */
                    if(reorientAtInteraction) {
                        cosines(theta, k);
                        for(d = 1; d <= Dims; d++) vel[k][d] = speed*theta[k][d]; 
                    }
                    /* Save connectivity matrix */
                    if(calcConn == 1) {
                        if(t > 0) {
                            if (firingTime >0) fprintf(out2,"\n");
                            firingTime += t;
                            fprintf(out2, "%f\t", firingTime);
                        }
                        fprintf(out2, "%d\t%d\t", p, k); 
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
    free_dmatrix(theta, 1, S, 1, Dims);
    free_dmatrix(pos, 1, S, 1, Dims);
    free_dmatrix(vel, 1, S, 1, Dims);
    free_ivector(Stime, 1, R);
    free_ivector(neigh, 1, S);

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

/* ----- Random Normal number generator ----- */
double ranNorm(void)
{
    double x, y;
    // Box-Muller method
    x = ranMT(); y = ranMT();
    return (double) sqrt(-2*log(x))*cos(2*PI*y);
}

/* ----- N-Dimensional directional cosines ----- */
void cosines(double **theta, int i)
{
    int d;
    double dum, sum = 0;

    for(d = 1; d <= Dims; d++) {
        dum =  ranNorm();
        sum += dum*dum; 
        theta[i][d] = dum;
    }
    for(d = 1; d <= Dims; d++) theta[i][d] /= sqrt(sum);
}

/* ------ Find all neighbors ----- */
void findNeighbors(double **pos, double **vel, int p, int *neigh)
{
    int i;
    /* Initialize neighborhood to zero */
    for(i = 1; i <= S; i++) neigh[i] = 0;
    /* Use appropriate function according to neighborhood */
    if(strcmp(neighborhood, "QNearest") == 0) qNearest(pos, p, neigh);
    else if(strcmp(neighborhood, "ConeOut") == 0) cone(pos, vel, p, neigh, "Out");
    else if(strcmp(neighborhood, "ConeIn") == 0) cone(pos, vel, p, neigh, "In");
}

/* ----- Find Q nearest neighbors ----- */
void qNearest(double **pos, int p, int *neigh)
{
    int i, j, foo;
    double dum = 2*L*L, dum1;
    double dis[9];
    double *dis2all;
    dis2all = dvector(1, S);
    dis2all[p] = dum;

    /* Calculate distances to all other oscillators */
    for(i = 1; i <= S; i++) {
        if(i == p) continue;
        dis2all[i] = edist(pos, p, i);
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
void cone(double **pos, double **vel, int p, int *neigh, char *direction)
{
    int j;
    double angle, dpos;
    for(j = 1; j <= S; j++) {
        if(j == p) continue;
        dpos = edist(pos, p, j);
        if(dpos <= r) {
            if (direction == "Out") angle = cangle(pos, vel, p, j);
            if (direction == "In") angle = cangle(pos, vel, j, p);
            if(angle <= alpha) neigh[j] = 1;            
        }
    }
}

/* ----- Euclidian distance between p & k ----- */
double edist(double **pos, int p, int k)
{
    int d;
    double dis = 0;
    for(d = 1; d <= Dims; d++) dis += (pos[p][d] - pos[k][d])*(pos[p][d] - pos[k][d]);
    dis = sqrt(dis); 
    return(dis);
}

/* ----- Angle between p->k and p's direction of motion ----- */
double cangle(double **pos, double **vel,int p,int k)
{
    int d;
    double angle = 0;

    /* Use scalar product between vel(p) and pos(k)-pos(p) */
    for(d = 1; d <= Dims; d++) angle += vel[p][d]*(pos[k][d] -pos[p][d]);
    angle /= (speed*edist(pos, p, k));
    angle = acos(angle); 
 
    return(angle);
}

/* ----- Calculate Order Parameter ----- */
double orderpar(double *phase)
{
    int i;
    double p = 0;
    for(i = 1; i <= S; i++) p += cos(2*PI*phase[i]);    
    p /= S;
    return(p);
}

/* ----- Movement in an unbounded environment (cyclic boundary conditions) ---- */
void unboundedMove(double **pos, double **vel, double t)
{
    int i, d;
    for(i = 1; i <= S; i++) {
        for(d = 1; d <= Dims; d++) {
            pos[i][d] += vel[i][d]*t;
            /* Periodical boundary conditions */
            while(pos[i][d] > L) pos[i][d] -= L;
            while(pos[i][d] < 0) pos[i][d] += L;
        }
    }
}

/* ----- Movement in a bounded environment (random reflections at wall) ---- */
void boundedMove(double **pos, double **vel, double **theta, double t)
{
    int i, d, q, bounceDim, dir;
    double tb, dum;

    while(t > 0) {
        tb = 2;
        /* Compute time for next boundary hit */
        for(i = 1; i <= S; i++) {
            for(d = 1; d <= Dims; d++) {
                if(vel[i][d] > 0) dum = (L - pos[i][d])/vel[i][d];
                if(vel[i][d] < 0) dum = (-pos[i][d])/vel[i][d];
                if(vel[i][d] == 0) dum = 2;
                if(dum < tb) {
                    tb = dum; q = i; bounceDim = d;
                }
            }
        }
        
        /* Boundary hit before next firing? */
        if(tb > t) tb = t;
        t -= tb;
        
        /* Update positions */
        for(i = 1; i <= S; i++) {
            for(d = 1; d <= Dims; d++) {
                pos[i][d] += vel[i][d]*tb;
            }
        }
        
        /* If wall is hit do random reflection */
        /* Distinguish close end (pos[q][bounceDim]=0) and far end (pos[q][bounceDim]=L) */
        dir = (vel[q][bounceDim] > 0) ? -1 : 1;
        if(t > 0) {
            cosines(theta, q);
            for(d = 1; d <= Dims; d++) {
                vel[q][d] = speed*theta[q][d];
                if(d == bounceDim) vel[q][d] = fabs(vel[q][d])*dir;
            }
        }
    }
}
