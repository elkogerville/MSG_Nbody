# MSG_Nbody

**MSG_Nbody** is a Python package offers an efficient, fully vectorized 3-dimensional NumPy implementation of the particle-particle N-body simulation algorithm by integrating the motion of stellar particles under their combined gravitational attraction. Initial conditions of different galaxy models in equilibrium are provided, including a Hernquist spherical galaxy and a simple disk galaxy. The algorithm for generating spherical galaxy initial conditions of different masses and scale lengths is also provided for further customizations. Yet, any set of initial conditions can be used as inputs to the simulation code, which will integrate their motions and save snapshot files directly to a directory. On a reasonably powerful personal computer, the code can support up to ~20,000 - 30,000 particles with runtimes on the order of a couple of days, using the numba compiler. Lowering the number of particles (N<15,000) will yield computation times of a couple minutes to a couple of hours. Therefore, this package aims to provide an accessible N-body simulation code in Python that is simple to set up and modify yet still simulates the effects of gravity with reasonable accuracy. The package also comes with a fully integrated Python toolkit to analyze simulation snapshots.

## Installation

You can install MSG_Nbody via pip:

pip install MSG-Nbody

The full documentation can be found on github at https://github.com/elkogerville/MSG_Nbody
