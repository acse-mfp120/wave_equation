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
Intialised with single radial disturbance         |  Initialised with two radial disturbances
:-------------------------:|:-------------------------:
![](https://github.com/acse-mfp120/wave_equation/blob/assets/visualised_results/4proc_output_360.0_360.0_40s.gif)  |  ![](https://github.com/acse-mfp120/wave_equation/blob/assets/visualised_results/4proc_dir_output_360.0_360.0_40s.gif)
