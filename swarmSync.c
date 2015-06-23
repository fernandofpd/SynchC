/* ---------- swarmSync.c --------- */
/* Swarm Synchronization of "numAgents" particles moving with speed "speed" in a square environment of length "length" */
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
double ranMT(void);
double ranNorm(void);
void cosines(double **theta, int i);
void findNeighbors(double **pos, double **vel, int p, int *neigh);
void qNearest(double **pos, int p, int *neigh);
void cone(double **pos, double **vel, int p, int *neigh, char *direction);
double edist(double **pos, int p, int k, int withShifts, double *shift);
double cangle(double **pos, double **vel, int p, int k, double dist, double *shift, int shiftSign);
double orderpar(double *phase);
void unboundedMove(double **pos, double **vel, double t);
void boundedMove(double **pos, double **vel, double **theta, double t);
void boundaryShifts();

/* ----------  Parameters ----------- */
/* Agent parameters */
double tau = 1.;                        // Oscillator period
double eps = 0.1;                       // Interaction strength, epsilon
int numAgents = 20;                     // Number of agents
double speed = 0.1;                     // Speed of the agents
/* Scenario parameters */
int dims = 2;                           // Number of Dimensions of the environment
double length = 400;                    // Length of environment (size is length^D)
int boundary = 0;                       // Bounded (=1) or Unbounded (=0) environment
/* Interaction parameters */
char neighborhood[128] = "QNearest";    // Neighborhood type (QNearest, ConeIn or ConeOut)
double alpha = 18;                      // Angle of Interaction (used in ConeIn & ConeOut)
double r = 400;                         // Raange of Interaction (used in ConeIn & ConeOut)
int Q = 1;                              // Number of Nearest Neighbors (used in QNearest)
double prob = 1;
/* Synchronization parameters */
double smin = 1.e-6;                    // Synchronization minimum, separation from perfect syncrhony
double Tmax = 1e7;                      // Maximum cycle count (censoring value)

/* ---------- Other Declarations ---------- */
int numShifts;
double **shifts;


/* ---------- Main ---------- */
int main(int argc,char **argv)
{
    /* ----- Variable Declarations and Inicializations -----*/
    int i, j, k, d, seed1 = 7, p;
    uint32 seed = 5;
    int calculateConnectivity, reorientAtInteraction, reorientAtFiring;    // Flags
    double t, firingTime, orderParam;
    double *phase, **theta, **pos,  **vel;                                 // Oscillators phases, motion angles, positions and velocities
    int *neigh;                                                            // Neighbor indices at each firing
    int numRuns = 1, *syncTimes, time;                                     // Number of Runs and syncronization times on each run
    FILE *out1, *out2;                                                     // Output files
    char filename[128] = "dat", filename2[128] = "conn";

    /* ----- Input parameters ----- */
    calculateConnectivity = opt_flag(&argc, argv, 2, "-c", "--conn");                  // Calculate connectivity (if flag is on)
    boundary = opt_flag(&argc, argv, 2, "-b", "--bounded");                            // Bounded environment (bounce at walls) (if flag is on) or has periodical boundary conditions
    reorientAtInteraction = opt_flag(&argc, argv, 1, "--reorientAtInteraction");       // Reorient upon receiving an interaction 
    reorientAtFiring = opt_flag(&argc, argv, 1, "--reorientAtFiring");                 // Reorient upon firing
    opt_int(&argc, argv, &Q, 2, "-Q", "--qnearest");                                   // Number of Nearest Neighbors
    opt_int(&argc, argv, &numAgents, 2, "-N", "--agents");                             // Total Number of Agents
    opt_int(&argc, argv, &numRuns, 2, "-R", "--runs");                                 // Number of Runs
    opt_int(&argc, argv, &dims, 2, "-D", "--dimensions");                              // Number of Runs    
    opt_int(&argc, argv, &seed1, 2, "-s", "--seed");                                   // Random Number Seed    
    opt_double(&argc, argv, &speed, 2, "-V", "--speed");                               // Speed of agents
    opt_double(&argc, argv, &eps, 2, "-e", "--epsilon");                               // Interaction parameter Epsilon
    opt_double(&argc, argv, &length, 2, "-length", "--length");                        // Length of the environment
    opt_double(&argc, argv, &alpha, 2, "-a", "--alpha");                               // Angle of Interaction (In/Out)
    opt_double(&argc, argv, &r, 2, "-r", "--radius");                                  // Radius of Interaction (In/Out)
    opt_double(&argc, argv, &Tmax, 2, "-T", "--tmax");                                 // Censoring threshold (stop run if system not synchronized before Tmax)
    opt_string(&argc, argv, filename, 2, "-f", "--filename");                          // Output filename
    opt_string(&argc, argv, neighborhood, 2, "-n", "--neighborhood");                  // Define the neighborhood model: {QNearest,ConeOut,ConeIn}

    /* ----- Initializations ----- */
    alpha = PI/180*alpha;

    out1 = fopen(filename, "w"); fclose(out1); out1 = fopen(filename, "a");
 
    phase = dvector(1, numAgents);
    theta = dmatrix(1, numAgents, 1, dims);
    pos = dmatrix(1, numAgents, 1, dims);
    vel = dmatrix(1, numAgents, 1, dims);
    syncTimes = ivector(1, numRuns);
    neigh = ivector(1, numAgents);

    boundaryShifts();

    seed = (uint32)(seed1); seedMT(seed);

    /* Loop over the Runs */
    for(j = 1; j <= numRuns; j++) {
        if(calculateConnectivity) {
            sprintf(filename2, "%sconn%d", filename, j);
            out2 = fopen(filename2, "w"); fclose(out2); out2 = fopen(filename2, "a");
        }

        orderParam = 0; time = 0; firingTime = 0;  // Reset

        /* Random initial values */
        for(i = 1; i <= numAgents; i++) {
            phase[i] = ranMT(); 
            cosines(theta, i);
            for(d = 1; d <= dims; d++) {
                pos[i][d] = length*ranMT();
                vel[i][d] = speed*theta[i][d];
            } 
        }

        /***** While not yet synchronized or before Tmax cycles have elapsed *****/
        while((orderParam < (1 -smin)) && time < Tmax) {
            /***** Update phases and perform movement *****/
            maxfind(&p, phase); t = 1 - phase[p];           // Index of next firing unit, p; time until next firing, t
            for(i = 1; i <= numAgents; i++) phase[i] += t;  // Update phases
            // Update positions until next firing
            if(t > 0) {
                if(boundary) boundedMove(pos, vel, theta, t);
                else unboundedMove(pos, vel, t);
            }

            /***** Firings and Interactions with Neighbors *****/
            /* When reference unit 1 fires calculate order parameter */
            if(p == 1) {
                orderParam = orderpar(phase);
                time++;
            }            
            phase[p] = 0;   // Reset phase
            /* If flag is active reorient upon firing */
            if(reorientAtFiring) {
                cosines(theta, p);
                for(d = 1; d <= dims; d++) vel[p][d] = speed*theta[p][d]; 
            }
            /* Find neighbors and update their phase */
            findNeighbors(pos, vel, p, neigh);
            for(k = 1; k <= numAgents; k++) {
                if(neigh[k] == 1) {  
                    phase[k] *= (1 + eps);
                    if(phase[k] >= 1) phase[k] = 1;
                    /* If flag is active reorient upon receiving an interaction */
                    if(reorientAtInteraction) {
                        cosines(theta, k);
                        for(d = 1; d <= dims; d++) vel[k][d] = speed*theta[k][d]; 
                    }
                    /* Save connectivity matrix */
                    if(calculateConnectivity == 1) {
                        if(t > 0) {
                            if (firingTime > 0) fprintf(out2, "\n");
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
        /* Save the synchronization times */
        syncTimes[j] = time;
        fprintf(out1, "%d\n", syncTimes[j]);

        if(calculateConnectivity) fclose(out2);
    }
 
    /* ----- Free memory ----- */
    fclose(out1);
    free_dvector(phase, 1, numAgents);
    free_dmatrix(theta, 1, numAgents, 1, dims);
    free_dmatrix(pos, 1, numAgents, 1, dims);
    free_dmatrix(vel, 1, numAgents, 1, dims);
    free_dmatrix(shifts, 1, numShifts, 1, dims); 
    free_ivector(syncTimes, 1, numRuns);
    free_ivector(neigh, 1, numAgents);

    return 0;
}

/* ----- Find index p corresponding to maximum phase ----- */
void maxfind(int *p, double *phase)
{
    int i;
    *p = 1;
    for(i = 2; i <= numAgents; i++) if(phase[i] > phase[*p]) *p = i;
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

    for(d = 1; d <= dims; d++) {
        dum =  ranNorm();
        sum += dum*dum; 
        theta[i][d] = dum;
    }
    for(d = 1; d <= dims; d++) theta[i][d] /= sqrt(sum);
}

/* ------ Find all neighbors ----- */
void findNeighbors(double **pos, double **vel, int p, int *neigh)
{
    int i;
    /* Initialize neighborhood to zero */
    for(i = 1; i <= numAgents; i++) neigh[i] = 0;
    /* Use appropriate function according to neighborhood */
    if(strcmp(neighborhood, "QNearest") == 0) qNearest(pos, p, neigh);
    else if(strcmp(neighborhood, "ConeOut") == 0) cone(pos, vel, p, neigh, "Out");
    else if(strcmp(neighborhood, "ConeIn") == 0) cone(pos, vel, p, neigh, "In");
}

/* ----- Find Q nearest neighbors ----- */
void qNearest(double **pos, int p, int *neigh)
{
    int i, j, foo;
    double dum = dims*pow(length, dims), dum1;
    double dis[9];
    double *dis2all;
    dis2all = dvector(1, numAgents);
    dis2all[p] = dum;

    /* Calculate distances to all other oscillators */
    for(i = 1; i <= numAgents; i++) {
        if(i == p) continue;
        dis2all[i] = edist(pos, p, i, 0, shifts[1]);
    }

    /* Find Q nearest neighbors */
    for(i = 1; i <= Q; i++) {
        dum1 = dum;
        for(j = 1; j <= numAgents; j++) {
            if(dis2all[j] < dum1) {
                foo = j;
                dum1 = dis2all[j];
            }
        }
        dis2all[foo] = dum;
        neigh[foo] = 1;
    }
    free_dvector(dis2all, 1, numAgents);
}

/* ----- Find Neighbors in a Cone ---- */
/* If direction == "Out" : p's neighbors will be the agents that lie inside its cone -- Cone of Influence */
/* If direction == "In" : p's neighbors will be the agents that see p inside their cone -- Cone of Vision */
void cone(double **pos, double **vel, int p, int *neigh, char *direction)
{
    int j, k;
    double angle, dist;
    for(j = 1; j <= numAgents; j++) {
        if(j == p) continue;
        for (k = 1; k <= numShifts; k++) {
            dist = edist(pos, p, j, 1, shifts[k]);
            if(dist <= r) {
                if (direction == "Out") angle = cangle(pos, vel, p, j, dist, shifts[k], 1);
                if (direction == "In") angle = cangle(pos, vel, j, p, dist, shifts[k], -1);
                if(angle <= alpha) {neigh[j] = 1; break;}
            }
        }
    }
}

/* ----- Euclidian distance between p & k ----- */
/* Shifts in the position of k are accepted through the array 'shift' (size 1xDim) if 'withShifts == 1' */
/* 'withShifts = 0' allows simple calculation of the shortest distance with periodic boundary conditions */
double edist(double **pos, int p, int k, int withShifts, double *shift)
{
    int d;
    double dis = 0, foo;

    for(d = 1; d <= dims; d++) {
        foo = abs(pos[k][d] - pos[p][d] + withShifts*shift[d]);
        // If periodical boundary conditions calculate shortest distance
        if (withShifts == 0 && boundary == 0 && foo > length/2) foo = length - foo;
        foo *= foo;
        dis += foo;
    }
    dis = sqrt(dis);
    return(dis);
}

/* ----- Angle between p->k and p's direction of motion ----- */
double cangle(double **pos, double **vel, int p, int k, double dist, double *shift, int shiftSign)
{
    int d;
    double angle = 0;

    /* Use scalar product between vel(p) and pos(k)-pos(p) */
    for(d = 1; d <= dims; d++) angle += vel[p][d]*(pos[k][d] - pos[p][d] + shiftSign*shift[d]);
    angle /= speed*dist;
    angle = acos(angle); 
 
    return(angle);
}

/* N-tuple of position shifts if periodic boundary conditions */
void boundaryShifts()
{
    int j, k, setSize = 3;
    double set[3] = {0, length, -length}, dum[dims];
    int i[dims];

    if (boundary == 0) numShifts = (int) pow(setSize, dims);
    else numShifts = 1;
    shifts = dmatrix(1, numShifts, 1, dims);
    
    // Initialize
    for(j = 0; j < dims; j++) {
        i[j] = 0;
        dum[j] = set[0];
        shifts[1][j+1] = set[0];
    }
    
    // Calculate all n-tuples if periodic boundary conditions
    k = 1;
    while(1) {
        if(k == numShifts) break;
        k++;
        for(j = 0; j < dims; j++) {
            i[j]++; 
            if(i[j] < setSize) { 
                dum[j] = set[i[j]];
                break;
            }
            i[j] = 0;
            dum[j] = set[0];
        }
        for(j = 0; j < dims; j++) shifts[k][j+1] = dum[j];
    }
}

/* ----- Calculate Order Parameter ----- */
double orderpar(double *phase)
{
    int i;
    double p = 0;
    for(i = 1; i <= numAgents; i++) p += cos(2*PI*phase[i]);    
    p /= numAgents;
    return(p);
}

/* ----- Movement in an unbounded environment (cyclic boundary conditions) ---- */
void unboundedMove(double **pos, double **vel, double t)
{
    int i, d;
    for(i = 1; i <= numAgents; i++) {
        for(d = 1; d <= dims; d++) {
            pos[i][d] += vel[i][d]*t;
            /* Periodical boundary conditions */
            while(pos[i][d] > length) pos[i][d] -= length;
            while(pos[i][d] < 0) pos[i][d] += length;
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
        for(i = 1; i <= numAgents; i++) {
            for(d = 1; d <= dims; d++) {
                if(vel[i][d] > 0) dum = (length - pos[i][d])/vel[i][d];
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
        for(i = 1; i <= numAgents; i++) {
            for(d = 1; d <= dims; d++) {
                pos[i][d] += vel[i][d]*tb;
            }
        }
        
        /* If wall is hit do random reflection */
        /* Distinguish close end (pos[q][bounceDim]=0) and far end (pos[q][bounceDim]=length) */
        dir = (vel[q][bounceDim] > 0) ? -1 : 1;
        if(t > 0) {
            cosines(theta, q);
            for(d = 1; d <= dims; d++) {
                vel[q][d] = speed*theta[q][d];
                if(d == bounceDim) vel[q][d] = fabs(vel[q][d])*dir;
            }
        }
    }
}
