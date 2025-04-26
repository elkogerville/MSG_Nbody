# MSG_Nbody

**MSG_Nbody** offers an efficient, fully vectorized and parallelized 3D NumPy implementation of the particle-particle N-body simulation algorithm. The simulation integrates the motion of stellar particles under their combined gravitational attraction in discretized timesteps. The acceleration computation is batch processed and compiled with Numba for fast calculation times. On a reasonably powerful personal computer, the code can support up to ~100,000 - 200,000 particles with runtimes on the order of a couple of days. Lowering the number of particles (N<60,000) will yield computation times from a couple minutes to a couple of hours. This package aims to provide an accessible N-body simulation code in Python that is simple to set up and modify yet still simulates the effects of gravity with reasonable accuracy. Additionally, initial conditions for different galaxy models in equilibrium are provided, including a spherical Hernquist galaxy and a simple disk galaxy. The algorithm for generating spherical galaxy initial conditions of different masses and scale lengths is also provided for further customizations, however, any set of initial conditions can be used as inputs to the simulation code. The package also comes with a fully integrated Python plotting library to analyze simulation snapshots.

## Installation

You can install MSG_Nbody via pip:

pip install MSG-Nbody

The full documentation can be found on github at https://github.com/elkogerville/MSG_Nbody
