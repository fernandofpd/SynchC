/* Fernando Perez-Diaz, February 2016           */

#include <math.h>

#include "tools.h"
#include "cokus.h" 
#include "variables.h"

/* ----- Find index p corresponding to maximum of vector (of length numAgents) ----- */
int maxFind(double *vector)
{
    int i, p = 1;
    for (i = 2; i <= numAgents; i++) if (vector[i] > vector[p]) p = i;
    return p;
}

/* ----- Find index p corresponding to minimum of vector (of length numAgents) ----- */
int minFind(double *vector)
{
    int i, p = 1;
    for (i = 2; i <= numAgents; i++) if (vector[i] < vector[p]) p = i;
    return p;
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

/* ----- N-Dimensional directional cosines to produce random orientations in N-D ----- */
void cosines(double *theta)
{
    int d;
    double dum, sum = 0;

    for (d = 1; d <= dims; d++) {
        dum =  ranNorm();
        sum += dum*dum; 
        theta[d] = dum;
    }
    for (d = 1; d <= dims; d++) theta[d] /= sqrt(sum);
}
