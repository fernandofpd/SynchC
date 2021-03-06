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
#include <string.h>

#include "neighborhood.h"
#include "movement.h"
#include "variables.h"
#include "tools.h"
#include "inputs.h"
#include "outputs.h"
#include "nrutil.h"
#include "cokus.h"

/* ----------- Function Declarations ----------- */
double orderpar(double *phase);
void initialize(double **phase_0, double *phase, double **tau_0, double *tau, double **pos_0, double **pos, double **vel_0, double **vel);
double phaseResponse(double phase);
int findNext(double *phase, double *tau);


/* ---------- Main ---------- */
int main(int argc, char **argv)
{
    /* ----- Variable Declarations and Initializations -----*/
    int i, j, k, p;
    uint32 seed = 5;
    int OUT_CONN, OUT_ORDPAR, OUT_INTERSPIKE;                                     // Flags related to output
    double dt, time, prevFiring, orderParam, lastTimePrint = 0;                   // Time to next firing, actual time, previous firing time of reference oscillator and order parameter
    double *phase, *tau, **pos,  **vel, **phase_0, **tau_0, **pos_0, **vel_0;     // Oscillators phases, natural period, positions and velocities
    int *neigh;                                                                   // Neighbor indices at each firing
    int Tcycles, Tsync;                                                           // Elapsed cycles and syncronization time
    FILE *tsyncFILE, *connFILE, *ordparFILE, *intspkFILE;                         // Output files
    char fphase[56] = "phases.txt", ftau[56] = "taus.txt", fpos[56] = "pos.txt", fvel[56] = "vel.txt";   // Input filenames
    // Output filenames
    char ftsync[56] = "dat";

    /* ----- Initializations ----- */
    inputOptions(argc, argv, &OUT_CONN, &OUT_ORDPAR, &OUT_INTERSPIKE, ftsync, fphase, ftau, fpos, fvel);
    boundaryShifts();
    seed = (uint32)(seed1); seedMT(seed);

    alpha = PI/180*alpha;

    tsyncFILE = fopen(ftsync, "w"); fclose(tsyncFILE);
 
    phase = dvector(1, numAgents);
    tau = dvector(1, numAgents);
    pos = dmatrix(1, numAgents, 1, dims);
    vel = dmatrix(1, numAgents, 1, dims);
    neigh = ivector(1, numAgents);

    phase_0 = dmatrix(1, numAgents, 1, 1);
    tau_0 = dmatrix(1, numAgents, 1, 1); 
    pos_0 = dmatrix(1, numAgents, 1, dims);
    vel_0 = dmatrix(1, numAgents, 1, dims);
    
    if (readInputFile(fphase, phase_0, numAgents, 1, 0, 1)) return EXIT_FAILURE;
    if (readInputFile(ftau, tau_0, numAgents, 1, 0, INFINITY)) return EXIT_FAILURE;
    if (readInputFile(fpos, pos_0, numAgents, dims, 0, length)) return EXIT_FAILURE;
    if (readInputFile(fvel, vel_0, numAgents, dims, 0, speed)) return EXIT_FAILURE;


    /* Loop over the Runs */
    for (j = 1; j <= numRuns; j++) {
        tsyncFILE = fopen(ftsync, "a");
        openFile(OUT_CONN, "conn", ftsync, &connFILE, j);
        openFile(OUT_ORDPAR, "ordPar", ftsync, &ordparFILE, j);
        openFile(OUT_INTERSPIKE, "interspike", ftsync, &intspkFILE, j);

	orderParam = 0; Tcycles = 0; Tsync = FLAG; time = 0; prevFiring = 0;     // Reset
        initialize(phase_0, phase, tau_0, tau, pos_0, pos, vel_0, vel);      // Initialize phases, positions and velocities


        /***** While not yet synchronized or before Tmax cycles have elapsed *****/
        while ((orderParam < (1 -smin) || STOPTMAX*(Tcycles < Tmax)) && (Tcycles < Tmax || STOPSYNC*(orderParam < (1 -smin)))) {
            p = findNext(phase, tau);                            // Index of next firing unit, p; 
            dt = (1 - phase[p])*tau[p];                          // Time until next firing, dt
            for (i = 1; i <= numAgents; i++)
                phase[i] += dt/tau[i];                           // Update phases of all units
            phase[p] = 0;                                        // Reset phase of firing unit

            if (dt > 0) {
                time += dt;                                      // Calculate absolute time of firing
 
                if (bounded) boundedMove(pos, vel, dt);          // Update positions until next firing
                else unboundedMove(pos, vel, dt);
            }

            if (p == 1) {
                orderParam = orderpar(phase);                    // When reference unit 1 fires calculate order parameter
                Tcycles++;
                if (OUT_ORDPAR) fprintf(ordparFILE, "%f\n", orderParam);    // Output Order Parameters
                if (OUT_INTERSPIKE && prevFiring != 0) fprintf(intspkFILE, "%f\n", time - prevFiring);  // Output Interspike Interval 
                prevFiring = time;
            }            

            if (REORIENT_FIRING) reorient(vel, p);               // If flag is active reorient upon firing

            findNeighbors(pos, vel, p, neigh);

            for (k = 1; k <= numAgents; k++) {
                if (neigh[k] == 1) {  
                    if (OUT_CONN) outConn(&connFILE, time, &lastTimePrint, p, k, phase); // Output Connectivity
                    if (REORIENT_INTERACTION) reorient(vel, k);  // If flag is active reorient upon receiving an interaction
                    phase[k] += phaseResponse(phase[k]);         // Update phases after interaction
                }
            }
            if (orderParam >= (1 -smin) && Tsync == FLAG) Tsync = Tcycles;
        }

        fprintf(tsyncFILE, "%d\n", Tsync); fclose(tsyncFILE);    // Save the synchronization times

        if (OUT_CONN) fclose(connFILE);
        if (OUT_ORDPAR) fclose(ordparFILE);
        if (OUT_INTERSPIKE) fclose(intspkFILE);
    }
 
    /* ----- Free memory ----- */
    free_dvector(phase, 1, numAgents);
    free_dvector(tau, 1, numAgents);
    free_dmatrix(pos, 1, numAgents, 1, dims);
    free_dmatrix(vel, 1, numAgents, 1, dims);
    free_dmatrix(shifts, 1, numShifts, 1, dims); 
    free_ivector(neigh, 1, numAgents);
    free_dmatrix(phase_0, 1, numAgents, 1, 1);
    free_dmatrix(tau_0, 1, numAgents, 1, 1);
    free_dmatrix(pos_0, 1, numAgents, 1, dims);
    free_dmatrix(vel_0, 1, numAgents, 1, dims);

    return 0;
}

/* ----- Calculate Order Parameter ----- */
double orderpar(double *phase)
{
    int i;
    double real = 0, im = 0;
    for (i = 1; i <= numAgents; i++) {
        real += cos(2*PI*phase[i]);
        im += sin(2*PI*phase[i]);
    }
    real /= numAgents;
    im /= numAgents;

    if (!strcmp(ordparFunc, "exp")) return sqrt(real*real + im*im);
    else return real;
}

/* ----- Update the phases of oscillators after receiving interaction ---- */
double phaseResponse(double phase)
{
    double dphase = 0;

    if (phase < refrac) dphase = 0;                   // No update during refractory period

    else if (!strcmp(responseFunc, "multiplicative")) dphase = eps*phase;  // Multiplicative update with factor eps

    else if (!strcmp(responseFunc, "sawtooth")) {
        if (phase < 0.5) dphase = -phase;             // Delay/inhibition with negative slope for phase < 0.5
        else dphase = 1 - phase;                      // Advance/excitation with negative slope otherwise
        dphase *= kappa;                              // Slope weight 0 < kappa <= 1
    }

    else if (!strcmp(responseFunc, "sine")) dphase = -kappa*sin(2*PI*phase)/(2*PI); // Like "sawtooth" but with variable slope

    if (!strcmp(responseFunc,"mixed")) {
        if (phase < refrac) dphase = -kappa*phase;    // Inhitbitory (delay) for phase < refrac
        else dphase = eps*phase;                      // Multiplicative otherwise
    }

    if ((phase + dphase) >= 1) dphase = 1 - phase;    // Prevent from exceeding threshold (effective negative slope)

    return dphase;
}


/* ----- Find index of next firing unit ----*/
int findNext(double *phase, double *tau)
{
    int i;
    double *deltaT;    // Time to next firing
    deltaT = dvector(1, numAgents);

    for (i = 1; i <= numAgents; i++) deltaT[i] = (1 - phase[i])*tau[i];
    
    free_dvector(deltaT, 1, numAgents);
    return minFind(deltaT);
}

/* ---- Initialize phases, positions and velocities ---- */
/* Default is random if input file not given or if value flagged on it */
/* Default value of tau is tauDefault */
void initialize(double **phase_0, double *phase, double **tau_0, double *tau, double **pos_0, double **pos, double **vel_0, double **vel)
{
    int i, d;
    double *theta;
    theta = dvector(1, dims);

    for (i = 1; i <= numAgents; i++) {
        if (phase_0[i][1] >= 0) phase[i] = phase_0[i][1];
        else phase[i] = ranMT();

        if (tau_0[i][1] != FLAG) tau[i] = tau_0[i][1];
        else tau[i] = tauDefault;
 
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
