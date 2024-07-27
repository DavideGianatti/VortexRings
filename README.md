February 2022

## Vortex Rings Simulation
This repository contains a simple Python script that simulates the evolution of vortex rings using the Biot-Savart law. The script generates, updates, and visualizes the positions and velocities of vortex tubes over a specified number of time steps. An example of vortex rings is smoke rings, which can interact in non trivial ways. Here below there is a sketch of what few steps of a simulation would look like.


### Parameters
The main parameters of the simulation are defined at the beginning of the script:

- Nt: Number of time evolution steps.
- dt: Time discretization.
- Nv: Number of tubes.
- M: Number of discretization points per tube.
- strength: Vorticity flux constant for each tube.
- length: Length of each tube.
- sigma: Width of each tube.
- position: Positions of discretization points for each tube.
- velocity: Velocities of discretization points.
- df: DataFrame to store the simulation data.

You can set up different scenarios by changing the inputs of the function **initialize_positions** at line 143.
