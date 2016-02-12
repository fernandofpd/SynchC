/* Fernando Perez-Diaz, February 2016           */

#include "variables.h"

/* ----------  Parameters ----------- */

/* Agent parameters */
double tau = 1.;                        // Oscillator period
double eps = 0.1;                       // Interaction strength, epsilon
int numAgents = 20;                     // Number of agents
double speed = 0.1;                     // Speed of the agents
int reorientInteraction = 0;
int reorientFiring = 0;

/* Scenario parameters */
int dims = 2;                           // Number of Dimensions of the environment
double length = 400;                    // Length of environment (size is length^D)
int bounded = 0;                        // Bounded (=1) or Unbounded (=0) environment

/* Interaction parameters */
char neighborhood[128] = "QNearest";    // Neighborhood type (QNearest, ConeIn or ConeOut)
double alpha = 18;                      // Angle of Interaction (used in ConeIn & ConeOut)
double r = 400;                         // Raange of Interaction (used in ConeIn & ConeOut)
int Q = 1;                              // Number of Nearest Neighbors (used in QNearest)
double prob = 1;                        // Probability of Interacting with Found Neighbors

/* Simulation parameters */
double smin = 1.e-6;                    // Synchronization minimum, separation from perfect syncrhony
double Tmax = 1e7;                      // Maximum cycle count (censoring value)
int numRuns = 1;                        // Number of runs
int seed1 = 7;                          // Seed for the random number generator
