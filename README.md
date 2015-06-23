# Swarm Synchronization
Swarm Synchronization of *NumAgents* particles moving with speed *vel* in a *D-dimensional* hypercube of length *L*. 
Particles are linear integrate and fire oscillators with period *tau*. 
When an oscillator's  phase reaches 1 it is reset to 0 and interacts with its neighbors. 
This interaction is a multiplicative update of the neighbor's phase by a factor *1 + epsilon*. 

Several type of neighborhoods are supported:  
- **QNearest**: p's neighbors are the Q closest particles.  
- **ConeOut**: p's neighbors will be the agents that lie inside its cone, i.e. *Cone of Influence*.  
- **ConeIn**: p's neighbors will be the agents that see p inside their cone, i.e *Cone of Vision*.  

**Bounded** and **unbounded** (cyclic) environments are supported. 
Reorientation to a random direction opposite to the wall's normal occur at wall hit for bounded environments. 
Additionally two other random reorientations are supported: **ReorientAtInteraction** and **ReorientAtFiring**, which occur upon firing or receiving an interaction respectively.

The system is updated either upon a wall-hit or a firing event. Therefore it is a pure continuous time simulation. 

The **order parameter** (which measures the degree of synchronization) is calculated at every firing event of the reference oscillator. 
The simulation runs until the order parameter reaches a certain threshold (*1-smin*) or *T<sub>max</sub>* cylces have elapsed. 
In the first case the number of cycles until synchronization, *T<sub>sync</sub>*, is saved, otherwise, the censoring at *T<sub>max</sub>* is indicated with an output of -1.  

Additionally, the connectivity matrix at each firing time can be saved.

Fernando Perez-Diaz. Ruediger Zillmer.
