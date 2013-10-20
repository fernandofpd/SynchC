/* randomd.h - double precision */
/* Routinen zur Erzeugung von Zufallszahlen */
/* Verwendet Routinen aus Numerical Recipes */
/* Volker Ahlers, Oktober 1998 */

void init_random(long *idum);
double gaussrandom(double mean,double stddev,long *idum);
double unirandom(double left,double right,long *idum);
double ran1(long *idum);
double gasdev(long *idum);
