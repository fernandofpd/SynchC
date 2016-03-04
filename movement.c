/* Fernando Perez-Diaz, February 2016           */

#include <math.h>

#include "movement.h"
#include "variables.h"
#include "tools.h"
#include "nrutil.h"

/* N-tuple of position shifts if periodic boundary conditions */
void boundaryShifts()
{
    int j, k, setSize = 3;
    double set[3] = {0, length, -length}, dum[dims];
    int i[dims];

    if (bounded == 0) numShifts = (int) pow(setSize, dims);
    else numShifts = 1;
    shifts = dmatrix(1, numShifts, 1, dims);
    
    // Initialize
    for (j = 0; j < dims; j++) {
        i[j] = 0;
        dum[j] = set[0];
        shifts[1][j+1] = set[0];
    }
    
    // Calculate all n-tuples if periodic boundary conditions
    for (k = 2; k <= numShifts; k++) {
        for (j = 0; j < dims; j++) {
            i[j]++; 
            if (i[j] < setSize) { 
                dum[j] = set[i[j]];
                break;
            }
            i[j] = 0;
            dum[j] = set[0];
        }
        for (j = 0; j < dims; j++) shifts[k][j+1] = dum[j];
    }
}


/* ----- Movement in an unbounded environment (cyclic boundary conditions) ---- */
void unboundedMove(double **pos, double **vel, double dt)
{
    int i, d;
    for (i = 1; i <= numAgents; i++) {
        for (d = 1; d <= dims; d++) {
            pos[i][d] += vel[i][d]*dt;
            /* Periodical boundary conditions */
            while (pos[i][d] > length) pos[i][d] -= length;
            while (pos[i][d] < 0) pos[i][d] += length;
        }
    }
}

/* ----- Movement in a bounded environment (random reflections at wall) ---- */
void boundedMove(double **pos, double **vel, double dt)
{
    int i, d, q, bounceDim, dir;
    double tb, dum;

    while (dt > 0) {
        tb = 2*tau;
        /* Compute time for next boundary hit */
        for (i = 1; i <= numAgents; i++) {
            for (d = 1; d <= dims; d++) {
                if (vel[i][d] > 0) dum = (length - pos[i][d])/vel[i][d];
                if (vel[i][d] < 0) dum = (-pos[i][d])/vel[i][d];
                if (vel[i][d] == 0) dum = 2*tau;
                if (dum < tb) {
                    tb = dum; q = i; bounceDim = d;
                }
            }
        }
        
        /* Boundary hit before next firing? */
        if (tb > dt) tb = dt;
        dt -= tb;

        /* Update positions */
        for (i = 1; i <= numAgents; i++) {
            for (d = 1; d <= dims; d++) {
                pos[i][d] += vel[i][d]*tb;
            }
        }
        
        /* If wall is hit do random reflection */
        /* Distinguish close end (pos[q][bounceDim]=0) and far end (pos[q][bounceDim]=length) */
        if (dt > 0) {
            dir = (vel[q][bounceDim] > 0) ? -1 : 1;
            reorient(vel, q);
            vel[q][d] = fabs(vel[q][d])*dir;
        }
    }
}

/* ---- Random Reorientation in N-Dimensions ---- */
void reorient(double **vel, int idx)
{
    int d;
    double *theta;
    theta = dvector(1, dims);

    cosines(theta);
    for (d = 1; d <= dims; d++) vel[idx][d] = speed*theta[d];

    free_dvector(theta, 1, dims);
}
