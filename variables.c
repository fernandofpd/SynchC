/* Fernando Perez-Diaz, February 2016           */

#include "variables.h"

/* ----------  Parameters ----------- */

/* Agent parameters */
double tauDefault = 1.;                    // Oscillator period
double eps = 0.1;                          // Interaction strength, epsilon, for multiplicative PRC
double kappa = 1;                          // Coupling strenght for delay-advance PRC
char responseFunc[50] = "multiplicative";  // Phase Response function ("multiplicative", "sawtooth", "sine")
double refrac = 0;                         // Refractory period (as a fraction of the phase)
int numAgents = 20;                        // Number of agents
double speed = 0.1;                        // Speed of the agents
int REORIENT_INTERACTION = 0;              // Flag, reorient upon receiving an interaction
int REORIENT_FIRING = 0;                   // Flag, reorient upon firing

/* Scenario parameters */
int dims = 2;                              // Number of Dimensions of the environment
double length = 400;                       // Length of environment (size is length^D)
int bounded = 0;                           // Bounded (=1) or Unbounded (=0) environment

/* Interaction parameters */
char neighborhood[128] = "QNearest";       // Neighborhood type (QNearest, ConeIn or ConeOut)
double alpha = 18;                         // Angle of Interaction (used in ConeIn & ConeOut)
double r = 400;                            // Raange of Interaction (used in ConeIn & ConeOut)
int Q = 1;                                 // Number of Nearest Neighbors (used in QNearest)
double prob = 1;                           // Probability of Interacting with Found Neighbors

/* Simulation parameters */
double smin = 1.e-6;                       // Synchronization minimum, separation from perfect syncrhony
double Tmax = 1e7;                         // Maximum cycle count (censoring value)
char ordparFunc[10] = "cos";               // Order parameter function 
int numRuns = 1;                           // Number of runs
int seed1 = 7;                             // Seed for the random number generator
int STOPTMAX = 0;                          // Stop simulation only when exactly Tmax cycles have elapsed
int STOPSYNC = 0;                          // Stop simulation only when synchronized, independently of the time it takes
