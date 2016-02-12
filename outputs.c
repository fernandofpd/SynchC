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
void outConn(FILE **connFILE, int *PRINT_TIME, double firingTime, double dt, int p, int k, double *phase)
{
    if (*PRINT_TIME) {
        if (firingTime != dt) fprintf(*connFILE, "\n"); // If it's not first line
        fprintf(*connFILE, "%f\t", firingTime);
        *PRINT_TIME = 0;
    }
    fprintf(*connFILE, "%d\t%d\t%f\t", p, k, phase[k]); 
}
