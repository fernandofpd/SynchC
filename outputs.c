/* Fernando Perez-Diaz, February 2016           */

#include <stdio.h>

/* ---- Open Output File if flag is on ---- */
void openFile(int flag, char *name, char *fname, FILE **file, int idx)
{
    char newName[50];
    if (flag) {
        sprintf(newName, "%s%s%d", fname, name, idx);
        *file = fopen(newName, "w");
    }
}

/* Save connectivity matrix */
void outConn(FILE **connFILE, double firingTime, double *lastTimePrint, int p, int k, double *phase)
{
    if (firingTime != *lastTimePrint) {
        if (*lastTimePrint != 0) fprintf(*connFILE, "\n"); // If it's not first line
        fprintf(*connFILE, "%f\t", firingTime);
        *lastTimePrint = firingTime;
    }
    fprintf(*connFILE, "%d\t%d\t%f\t", p, k, phase[k]); 
}
