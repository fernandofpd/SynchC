/* Fernando Perez-Diaz, February 2016           */

#include <stdio.h>
#include <stdarg.h>

#include "inputs.h"
#include "variables.h"
#include "options.h"

/* ----- Read a matrix from an input file and save it into the corresponding 2D array ---- */
/* Checks that input is of size = numRows x numCols */
/* Each value must lie in the [minVal, maxVal] interval or be equal to FLAG (-1 = FLAG < minVal < maxVal) */
int readInputFile(char *inputFilename, double **matrix, int numRows, int numCols, double minVal, double maxVal)
{
    int col = 0, row = 1;
    char next;
    double dum;
    char ERROR[256];
    FILE *fileIN;

    sprintf(ERROR, "ERROR: %s must be a  %d x %d matrix with values in the range [%f, %f] or %d to flag random initialization", inputFilename, numRows, numCols, minVal, maxVal, FLAG);
    fileIN = fopen(inputFilename, "r");

    if (fileIN) {
        while (fileIN && !feof(fileIN)) {
            int count = fscanf(fileIN, "%lf%c", &dum, &next);
            if (count > 0) {
                if ((dum < minVal || dum > maxVal) && dum != FLAG) {printf("%s\n", ERROR); return 1;}
                col++;
                if (col <= numCols) matrix[row][col] = dum;
                else {printf("%s\n", ERROR); return 1;}
                if (count == 1 || next == '\n') {
                    if (row > numRows) {printf("%s\n", ERROR); return 1;}
                    if (col < numCols) {printf("%s\n", ERROR); return 1;}
                    row++; col = 0; 
                }
            }
        }
        if (row <= numRows) {printf("%s\n", ERROR); return 1;}
        fclose(fileIN);
    }
    else for (row = 1; row <= numRows; row++) for (col = 1; col <= numCols; col++) matrix[row][col] = FLAG;

    return 0;
}

/* ----- Input Options ----- */
void inputOptions(int argc, char *argv[], int *outConn, int *outOrdPar, int *outInterspike, char *ftsync, char  *fphase, char *ftau, char *fpos, char *fvel)
{
    /* Define Parameters */
    opt_int(&argc, argv, &numAgents, 2, "-N", "--agents");                       // Total Number of Agents
    opt_int(&argc, argv, &dims, 2, "-D", "--dimensions");                        // Number of Runs    
    opt_double(&argc, argv, &speed, 2, "-V", "--speed");                         // Speed of agents
    opt_double(&argc, argv, &eps, 2, "-e", "--epsilon");                         // Interaction parameter Epsilon (multiplicative PRC)
    opt_double(&argc, argv, &kappa, 2, "-k", "--kappa");                         // Interaction parameter Kappa (delay-response PRC)
    opt_double(&argc, argv, &refrac, 1, "--refractory");                         // Refractory period
    opt_double(&argc, argv, &length, 2, "-L", "--length");                       // Length of the environment
    opt_double(&argc, argv, &Tmax, 2, "-T", "--tmax");                           // Censoring threshold (stop run if system not synchronized before Tmax)
 
    opt_int(&argc, argv, &Q, 2, "-Q", "--qnearest");                             // Number of Nearest Neighbors
    opt_double(&argc, argv, &alpha, 2, "-a", "--alpha");                         // Angle of Interaction (In/Out)
    opt_double(&argc, argv, &r, 2, "-r", "--radius");                            // Radius of Interaction (In/Out)
    opt_double(&argc, argv, &r, 2, "-p", "--prob");                              // Probability of interaction with neighbors

    opt_int(&argc, argv, &numRuns, 2, "-R", "--runs");                           // Number of Runs
    opt_int(&argc, argv, &seed1, 2, "-s", "--seed");                             // Random Number Seed
    STOPSYNC = opt_flag(&argc, argv, 1, "--stopsync");                           // Stop simulation only when synchronized, independently of the time it takes
    STOPTMAX = opt_flag(&argc, argv, 1, "--stoptmax");                           // Stop simulation only when exactly Tmax cycles have elapsed

    /* Define characteristics of movement, interaction or model */
    bounded = opt_flag(&argc, argv, 2, "-b", "--bounded");                       // Bounded environment (bounce at walls) (if flag is on) or has periodical boundary conditions
    REORIENT_INTERACTION = opt_flag(&argc, argv, 1, "--REORIENT_INTERACTION");   // Reorient upon receiving an interaction 
    REORIENT_FIRING = opt_flag(&argc, argv, 1, "--REORIENT_FIRING");             // Reorient upon firing
    opt_string(&argc, argv, neighborhood, 2, "-n", "--neighborhood");            // Define the neighborhood model: {QNearest,ConeOut,ConeIn}
    opt_string(&argc, argv, responseFunc, 1, "--response");                      // Define the phase response function: {multiplicative,sawtooth,sine}
    opt_string(&argc, argv, ordparFunc, 1, "--ordparfun");                       // Define the order parameter function: {cos,exp}

    /* Set desired outputs */
    *outConn = opt_flag(&argc, argv, 2, "-c", "--conn");                         // Calculate connectivity (if flag is on)
    *outOrdPar = opt_flag(&argc, argv, 1, "--ordpar");                           // Output order parameter (if flag is on)
    *outInterspike = opt_flag(&argc, argv, 1, "--interspike");                   // Output firing times (if flag is on)

    /* Set input or output file names */
    opt_string(&argc, argv, ftsync, 2, "-f", "--filename");                      // Output filename
    opt_string(&argc, argv, fphase, 1, "--phases");                              // Name of input file with agents initial phases
    opt_string(&argc, argv, ftau, 1, "--taus");                                  // Name of input file with agents initial phases
    opt_string(&argc, argv, fpos, 1, "--positions");                             // Name of input file with agents initial positions
    opt_string(&argc, argv, fvel, 1, "--velocities");                            // Name of input file with agents initial velocities 
} 
