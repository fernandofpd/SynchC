/* Fernando Perez-Diaz, February 2016           */

void findNeighbors(double **pos, double **vel, int p, int *neigh);
void qNearest(double **pos, int p, int *neigh);
void cone(double **pos, double **vel, int p, int *neigh, char *direction);
double edist(double **pos, int p, int k, double *shift);
double cangle(double **pos, double **vel, int p, int k, double dist, double *shift, int shiftSign);
