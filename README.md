# Repulser GSGP

Implementation of GSGP using the concept of semantic repulsors. 

## Idea
Detect overfitting solutions during execution of Geometric Semantic Genetic Programming and classify these as repulsers. In further generations transform standard GSGP into a multi-objective system, with objectives made of the fitness of the solutions and a distance metric to the repulsers. By maximizing the distance to repulsers and still looking for improved fitness, the system should help to limit overfitting of long running executions.
* based on [GSGP](http://gsgp.sourceforge.net/)
* modified using NSGA-II and a dynamic number of repulsors, to which the distance has to be maximized.

### Execution/Compilation Notes

* use g++ -std=c++11 (needed for tuple type)
