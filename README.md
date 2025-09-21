
# DEM Ball Drop Simulation

A minimal, verified Discrete Element Method simulation in Fortran 90. It models a ball dropped from a height `h`, calculating its trajectory, velocities, forces, and energy throughout the fall.

## Key Features
- **Physics:** Models motion under gravity using Newton's second law.
- **Output:** Generates set of plots describing the movement of the ball.
- **Verification:** The results have been verified for correct dynamics and energy conservation.

## Quick Start
1.  **Compile** with `gfortran`:
    ```bash
    gfortran -o dem_simulation_1D.f90 -O2
    ```
2.  **Run** the executable:
    ```bash
    ./dem_simulation_1D
    ```
3.  **Analyze** the plot in the `Plots/` directory.

## Configuration
Modify parameters directly in the `main.f90` source file:
- `initial_height`
- `gravity`
- `timestep (dt)`
- `particle_mass`
