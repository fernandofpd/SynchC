/* Fernando Perez-Diaz, February 2016           */

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "neighborhood.h"
#include "variables.h"
#include "nrutil.h"
#include "tools.h"

/* ------ Find all neighbors ----- */
void findNeighbors(double **pos, double **vel, int p, int *neigh)
{
    int i;
    /* Initialize neighborhood to zero */
    for (i = 1; i <= numAgents; i++) neigh[i] = 0;
    /* Use appropriate function according to neighborhood */
    if (strcmp(neighborhood, "QNearest") == 0) qNearest(pos, p, neigh);
    else if (strcmp(neighborhood, "ConeOut") == 0) cone(pos, vel, p, neigh, "Out");
    else if (strcmp(neighborhood, "ConeIn") == 0) cone(pos, vel, p, neigh, "In");
    /* Neighbors only with certain probability */
    for (i = 1; i <= numAgents; i++) {
        if (neigh[i] == 1) neigh[i] = ranMT() <= prob;
    }
}

/* ----- Find Q nearest neighbors ----- */
void qNearest(double **pos, int p, int *neigh)
{
    int i, k, q;
    double max = dims*pow(length, dims), dist, dum = max;
    double *dis2all;
    dis2all = dvector(1, numAgents);
    dis2all[p] = max;

    /* Calculate distances to all other oscillators */
    for (i = 1; i <= numAgents; i++) {
        if (i == p) continue;
        dum = max;
        for (k = 1; k <= numShifts; k++) {
            dist = edist(pos, p, i, shifts[k]);
            if (dist < dum) dum = dist;
        }
        dis2all[i] = dum;
    }
    
    /* Find Q nearest neighbors */
    for (i = 1; i <= Q; i++) {
        q = minFind(dis2all);
        dis2all[q] = max;
        neigh[q] = 1;
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
    for (j = 1; j <= numAgents; j++) {
        if (j == p) continue;
        for (k = 1; k <= numShifts; k++) {
            dist = edist(pos, p, j, shifts[k]);
            if (dist <= r) {
                if (strcmp(direction, "Out") == 0) angle = cangle(pos, vel, p, j, dist, shifts[k], 1);
                if (strcmp(direction, "In") == 0) angle = cangle(pos, vel, j, p, dist, shifts[k], -1);
                if (angle <= alpha) {neigh[j] = 1; break;}
            }
        }
    }
}

/* ----- Euclidian distance between p & k ----- */
/* Shifts in the position of k are accepted through the array 'shift' (size 1xDim) */
double edist(double **pos, int p, int k, double *shift)
{
    int d;
    double dis = 0, foo;

    for (d = 1; d <= dims; d++) {
        foo = fabs(pos[k][d] - pos[p][d] + shift[d]);
        dis += foo*foo;
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
    for (d = 1; d <= dims; d++) angle += vel[p][d]*(pos[k][d] - pos[p][d] + shiftSign*shift[d]);
    angle /= speed*dist;
    angle = acos(angle); 
 
    return(angle);
}

