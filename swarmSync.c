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

#include "neighborhood.h"
#include "movement.h"
#include "variables.h"
#include "tools.h"
#include "inputs.h"
#include "outputs.h"
#include "nrutil.h"
#include "options.h"
#include "cokus.h"

/* ----------- Function Declarations ----------- */
double orderpar(double *phase);
void initialize(double **phase_0, double *phase, double **pos_0, double **pos, double **vel_0, double **vel);
double phaseResponse(double phase);


/* ---------- Main ---------- */
int main(int argc,char **argv)
{
    /* ----- Variable Declarations and Initializations -----*/
    int i, j, k, p, PRINT_TIME;
    uint32 seed = 5;
    int OUT_CONN, OUT_ORDPAR, OUT_INTERSPIKE;                                     // Flags related to output
    double dt, time, prevFiring, orderParam;                                      // Time to next firing, actual time, previous firing time of reference oscillator and order parameter
    double *phase, **pos,  **vel, **phase_0, **pos_0, **vel_0;                    // Oscillators phases, positions and velocities
    int *neigh;                                                                   // Neighbor indices at each firing
    int Tsync;                                                                    // Syncronization times on each run
    FILE *tsyncFILE, *connFILE, *ordparFILE, *intspkFILE;                         // Output files
    char fphase[56] = "phases.txt", fpos[56] = "pos.txt", fvel[56] = "vel.txt";   // Input filenames
    // Output filenames
    char ftsync[56] = "dat";

    /* ----- Initializations ----- */
    inputOptions(argc, argv, &OUT_CONN, &OUT_ORDPAR, &OUT_INTERSPIKE, ftsync, fphase, fpos, fvel);
    boundaryShifts();
    seed = (uint32)(seed1); seedMT(seed);

    alpha = PI/180*alpha;

    tsyncFILE = fopen(ftsync, "w");
 
    phase = dvector(1, numAgents);
    pos = dmatrix(1, numAgents, 1, dims);
    vel = dmatrix(1, numAgents, 1, dims);
    neigh = ivector(1, numAgents);

    phase_0 = dmatrix(1, numAgents, 1, 1);
    pos_0 = dmatrix(1, numAgents, 1, dims);
    vel_0 = dmatrix(1, numAgents, 1, dims);
    
    if (readInputFile(fphase, phase_0, numAgents, 1, 0, 1)) return EXIT_FAILURE;
    if (readInputFile(fpos, pos_0, numAgents, dims, 0, length)) return EXIT_FAILURE;
    if (readInputFile(fvel, vel_0, numAgents, dims, 0, speed)) return EXIT_FAILURE;


    /* Loop over the Runs */
    for (j = 1; j <= numRuns; j++) {
        openFile(OUT_CONN, "conn", ftsync, &connFILE, j);
        openFile(OUT_ORDPAR, "ordPar", ftsync, &ordparFILE, j);
        openFile(OUT_INTERSPIKE, "interspike", ftsync, &intspkFILE, j);

        orderParam = 0; Tsync = 0; time = 0; prevFiring = 0;     // Reset
        initialize(phase_0, phase, pos_0, pos, vel_0, vel);      // Initialize phases, positions and velocities


        /***** While not yet synchronized or before Tmax cycles have elapsed *****/
        while ((orderParam < (1 -smin)) && Tsync < Tmax) {

            p = maxFind(phase); dt = 1 - phase[p];               // Index of next firing unit, p; time until next firing, t
            for (i = 1; i <= numAgents; i++) phase[i] += dt;     // Update phases of other units
            phase[p] = 0;                                        // Reset phase of firing unit

            if (dt > 0) {
                time += dt;                                      // Calculate absolute time of firing
 
                if (bounded) boundedMove(pos, vel, dt);          // Update positions until next firing
                else unboundedMove(pos, vel, dt);
            }

            if (p == 1) {
                orderParam = orderpar(phase);                    // When reference unit 1 fires calculate order parameter
                Tsync++;
                if (OUT_ORDPAR) fprintf(ordparFILE, "%f\n", orderParam);    // Output Order Parameters
                if (OUT_INTERSPIKE && prevFiring != 0) fprintf(intspkFILE, "%f\n", time - prevFiring);  // Output Interspike Interval 
                prevFiring = time;
            }            

            if (REORIENT_FIRING) reorient(vel, p);               // If flag is active reorient upon firing

            findNeighbors(pos, vel, p, neigh);

            PRINT_TIME = dt > 0;                                 //flag to print out the firing time if OUT_CONN = 1
            for (k = 1; k <= numAgents; k++) {
                if (neigh[k] == 1) {  
                    if (OUT_CONN) outConn(&connFILE, &PRINT_TIME, time, dt, p, k, phase); // Output Connectivity
                    if (REORIENT_INTERACTION) reorient(vel, k);  // If flag is active reorient upon receiving an interaction
                    phase[k] += phaseResponse(phase[k]);         // Update phases after interaction
                }
            }
        }
        if (Tsync >= Tmax) Tsync = FLAG;                         // If not synchronized in Tmax cycles, censor data

        fprintf(tsyncFILE, "%d\n", Tsync);                       // Save the synchronization times

        if (OUT_CONN) fclose(connFILE);
        if (OUT_ORDPAR) fclose(ordparFILE);
        if (OUT_INTERSPIKE) fclose(intspkFILE);
    }
 
    /* ----- Free memory ----- */
    fclose(tsyncFILE);
    free_dvector(phase, 1, numAgents);
    free_dmatrix(pos, 1, numAgents, 1, dims);
    free_dmatrix(vel, 1, numAgents, 1, dims);
    free_dmatrix(shifts, 1, numShifts, 1, dims); 
    free_ivector(neigh, 1, numAgents);
    free_dmatrix(phase_0, 1, numAgents, 1, 1);
    free_dmatrix(pos_0, 1, numAgents, 1, dims);
    free_dmatrix(vel_0, 1, numAgents, 1, dims);

    return 0;
}

/* ----- Calculate Order Parameter ----- */
double orderpar(double *phase)
{
    int i;
    double p = 0;
    for (i = 1; i <= numAgents; i++) p += cos(2*PI*phase[i]);    
    p /= numAgents;
    return(p);
}

/* ------ Update the phases of oscillators after receiving interaction ---- */
double phaseResponse(double phase)
{
    double dphase = 0;

    if (phase < refrac) dphase = 0;    // No update during refractory period

    else if (!strcmp(responseFunc, "multiplicative")) dphase = eps*phase;

    else if (!strcmp(responseFunc, "sawtooth")) {
        if (phase < 0.5) dphase = -phase;
        else dphase = 1 - phase; 
    }

    else if (!strcmp(responseFunc, "sine")) dphase = -sin(2*PI*phase)/(2*PI);

    if ((phase + dphase) >= 1) dphase = 1 - phase;    // Prevent from exceeding threshold

    return dphase;
}

/* ---- Initialize phases, positions and velocities ---- */
/* Default is random if input file not given or if value flagged on it */
void initialize(double **phase_0, double *phase, double **pos_0, double **pos, double **vel_0, double **vel)
{
    int i, d;
    double *theta;
    theta = dvector(1, dims);

    for (i = 1; i <= numAgents; i++) {
        if (phase_0[i][1] >= 0) phase[i] = phase_0[i][1];
        else phase[i] = ranMT();
 
        cosines(theta);

        for (d = 1; d <= dims; d++) {
            if (pos_0[i][d] != FLAG) pos[i][d] = pos_0[i][d];
            else pos[i][d] = length*ranMT();

            if (vel_0[i][d] >= 0) vel[i][d] = vel_0[i][d];
            else vel[i][d] = speed*theta[d];
        }
    }
    free_dvector(theta, 1, dims);
}
