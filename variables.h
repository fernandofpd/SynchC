/* Fernando Perez-Diaz, February 2016           */

#define FLAG -1

/* ----------  Parameters ----------- */

/* Agent parameters */
extern double tau;               // Oscillator period
extern double eps;               // Interaction strength, epsilon
extern int numAgents;            // Number of agents
extern double speed;             // Speed of the agents
extern int reorientInteraction;
extern int reorientFiring;

/* Scenario parameters */
extern int dims;                 // Number of Dimensions of the environment
extern double length;            // Length of environment (size is length^D)
extern int bounded;              // Bounded (=1) or Unbounded (=0) environment

/* Interaction parameters */
extern char neighborhood[128];   // Neighborhood type (QNearest, ConeIn or ConeOut)
extern double alpha;             // Angle of Interaction (used in ConeIn & ConeOut)
extern double r;                 // Raange of Interaction (used in ConeIn & ConeOut)
extern int Q;                    // Number of Nearest Neighbors (used in QNearest)
extern double prob;              // Probability of Interacting with Found Neighbors

/* Simulation parameters */
extern double smin;              // Synchronization minimum, separation from perfect syncrhony
extern double Tmax;              // Maximum cycle count (censoring value)
extern int numRuns;              // Number of runs
extern int seed1;                // Seed for the random number generator

/* ---------- Other Declarations ---------- */
int numShifts;
double **shifts;
