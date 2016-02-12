/* Fernando Perez-Diaz, February 2016           */

void unboundedMove(double **pos, double **vel, double dt);
void boundedMove(double **pos, double **vel, double dt);
void boundaryShifts();
void reorient(double **vel, int idx);
