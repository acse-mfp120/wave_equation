# 2D Wave Equation Solver

Solving the wave equation:

$$ \frac{\partial^2u}{\partial t^2} = c^2\nabla^2u $$

using a finite difference discretisation.

A C++ implementation parallelised with MPI, run and tested on an HPC cluster. Features domain decomposition, non-blocking point to point and collective comms and custom MPI variables.

## Configurations:
Set config.txt to choose between:
- Boundary condition (Dirichlet and Neumann)
- Fixed or periodic boundaries
- Grid size and resolution
- Time resolution
- Wavespeed, $c$
- Initial, radial disturbances
- Max runtime

## Example Outputs



## Performance
