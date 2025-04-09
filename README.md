# Merging Simulator for Galaxies [MSG]
## MSG_Nbody


<p align="center">
  <img src="ANIMATIONS/ezgif-6475adeb9afc4b.gif" width="45%">
  <img src="ANIMATIONS/ezgif-65d94a97da5a3c.gif" width="45%">
</p>


<div style="clear: both; margin-top: 20px;">
  
#### N-body simulation code using the particle-particle algorithm for studying galaxy mergers with particles on the order of ~10<sup>4</sup>

## Installation
```
$ pip install MSG_Nbody
```

<div align="justify"> 

## see [Documentation](#documentation)
For a guide on how to setup a simulation from start to finish, see the [MSG_Nbody Documentation Notebook]

## Uses
The MSG_Nbody Python package offers an efficient, fully vectorized 3-dimensional NumPy implementation of the particle-particle N-body simulation algorithm which integrates the motion of stellar particles through space under their mutual gravitational attraction. Initial conditions of different galaxy models in equilibrium are provided, including a Hernquist spherical galaxy and a simple disk galaxy. The algorithm for generating spherical galaxy initial conditions of different masses and scale lengths is also provided for further customizations. Yet, any set of initial conditions can be used as inputs to the simulation code, which will integrate their motions and output snapshot files directly to a directory. On a reasonably powerful personal computer, the code can support up to ~20,000 - 30,000 particles with runtimes on the order of a couple of days, using the numba compiler. Lowering the number of particles (N<15,000) will yield computation times on the order of a couple minutes to a couple of hours. Therefore, the purpose of this package is to provide an accessible N-body simulation code in python that is simple to set up and modify yet still simulates the effects of gravity with reasonable accuracy. The package also comes with a fully integrated analysis toolkit to analyze simulation snapshots.

## Documentation and How to Use
For an in-depth documentation review please see the 'Documentation_and_Startup_Guide' jupyter notebook in the 'Tests' directory. This is an overview of how to set up and run an N-body simulation between two colliding disk galaxies with 6,000 total particles. In the 'Tests' folder there are numerous programs for the user to try, including jupyter notebook guides as well as Python programs to run from the terminal. The 'Generate_Spherical_Galaxy_Conditions' notebook demonstrates how to use the spherical_galaxy() function to generate initial conditions for a spherically symmetric galaxy. The N-body simulation code is flexible, however, and can take any set of initial conditions in the form of NumPy arrays as inputs. The code below demonstrates a simple 1:1 merger between two Hernquist spherical galaxies.
```python
from MSG_Nbody import *
# create spherical galaxy initial conditions, with a 2:1 mass ratio
galaxy_1 = spherical_galaxy(N=15000, mass=1, scale_length=1)
galaxy_2 = spherical_galaxy(N=13000, mass=1, scale_length=1)
# seperate into position, velocity, and mass arrays
gal1_pos, gal1_vel, gal1_mas = galaxy_1[:,0:3], galaxy_1[:,3:6], galaxy_1[:,6:7]
gal2_pos, gal2_vel, gal2_mas = galaxy_2[:,0:3], galaxy_2[:,3:6], galaxy_2[:,6:7]
# move perturber galaxy away from host galaxy and set it on collision path
gal2_pos += [40, 40, 0]
gal2_vel += [-0.2, -0.12, 0]
# append positions, velocities, and masses into 3 master arrays
pos_list = [gal1_pos, gal2_pos]
vel_list = [gal1_vel, gal2_vel]
mass_list = [gal1_mas, gal2_mas]
positions, velocities, masses = concatenate_initial_conditions(pos_list, vel_list, mass_list,
# run N-body simulation
dt = 0.01
timesteps = 5000
MSGnbody(positions, velocities, masses, dt, timesteps)
```
## The N-Body Problem


The N-body problem in astrophysics attempts to solve the motion of $N$ bodies through space under their mutual gravitational attraction. For a system with $N = 2$ bodies, there exists an analytical solution to their trajectories, allowing for accurate predictions of their positions and velocities at some future time $t$. The problem arises when $N \geq$ 3, where the chaotic nature of the system results in no solvable analytic solution. Chaotic systems are not random in nature but instead are characterized by having a unique solution to every set of initial conditions. These systems are highly sensitive to changes in initial conditions, where seemingly small fluctuations can lead to highly divergent solutions. For small numbers of $N$, constraints can be made to approximate the trajectories accurately. Yet when studying globular clusters or galaxies, $N \simeq 10^6-10^{11}$, therefore complicating calculations a great deal. Thus, the lack of an analytical solution requires a numerical approach for predicting the orbits of systems with large $N$. Numerical approximations integrate the equations of motion of each particle in discrete timesteps $\Delta t$, and then recursively use the previous set of positions and velocities to compute the next timestep. 



In this simulation, all particles are assumed to be collisionless, baryonic, stellar masses which are affected solely by gravitational forces. Collisions can be ignored completely when simulating stellar particles only due to the relaxation time of the system, which is defined as how long it takes for some star's trajectory to be significantly perturbed by the other stars in the system (equation 1). For a typical galaxy containing $N = 10^{11}$ stars and a crossing time of ~ $10^8$ years (average time for a star to cross the galaxy; $t_{cross} = \frac{R}{v}$), the relaxation time is orders of magnitude larger than a typical simulation timescale of a couple Gigayears. Thus, collisions can safely be removed as long as gas is not included in the simulation, which is collisional. The N-body code therefore assumes the only force operating on the particles is from their mutual gravitational attraction. Given a set of particles with initial positions and velocities, the next timestep is computed by brute force using the leap-frog algorithm. For each particle, the gravitational acceleration acting onto it must be calculated by summing up the individual particle-particle contributions from all the other stars. Thus, the gravitational acceleration $g_i$ onto a particle $p_{i}$ can be expressed as a sum over all the other particles $j$, where $r$ represents the particle's positional vector in 3D space and $\epsilon$, the softening length (equation 2). $\epsilon$ ensures the effects of close encounters are smoothed, and that dividing by zero does not occur. Its value is determined by the number of particles based on a relation derived by Dehnen (2001) for Plummer spheres (equation 3), and also serves as the simulation resolution. Close encounters between particles with a distance smaller than $\epsilon$ cannot be resolved.
<br>
$$t_{relax}\simeq\frac{N}{8lnN}\frac{R}{v}  \qquad (1) \qquad g_{i} = G\sum_{j}^{N}\frac{m_{j}[r_{j}-r_{i}]}{[|r_{j}-r_{i}|^2 + \epsilon ^2]^{3/2}} \qquad (2) \qquad \epsilon = 0.017 \left[ \frac{N}{10^5} \right]^{-0.23} \qquad (3)$$

Because this algorithm calculates the gravitational force from each particle onto each particle resulting in O($N^2$) calculations, the number of particles must be kept down. Thus the dark matter halo, central bulge and black hole, as well as gas particles, which are important components of galaxies are completely omitted in these simulations.

<figure>
  <img src="ANIMATIONS/sim2_pvd_t0_(0_-1_0).png" width="1000" align = 'center'>
</figure>


## N-Body Particle-Particle Algorithm
Once the gravitational acceleration onto each particle is computed for a given timestep using equation 2, the positions and velocities of the next timestep can then be calculated using the standard kinematic equations of motion (equations 4 and 5). The leap-frog algorithm computes the velocities and positions at interleaved timesteps where the velocities are calculated at half timesteps before and after computing the new positions. This creates a ’kick,’ ’drift,’ ’kick’ method conserving Energy to the second order and is a good trade-off between accuracy and computational efficiency. The new positions are then used to calculate a new set of accelerations, continuing the cycle endlessly.
<br>
$$v_{t+\frac{1}{2}} = v_{t} + g_{t} \frac{\Delta t}{2} \qquad (4) \qquad \qquad x_{t+1} = x_{t} + v_{t+\frac{1}{2}} \Delta t \qquad (5) \qquad \qquad \phi_{i} = \sum_{j} \frac{Gm_{j}}{|r_{j}-r_{i} + \epsilon|} \qquad (6)$$
<br>
The integrator saves the phase space coordinates $x,y,z,v_{x},v_{y},v_{z}$, and normalized potential $\phi_{i}$ (equation 6) of each particle every 10 timesteps as a $Nx7$ matrix. Moreover, in all simulations runs model units are assumed, where the gravitational constant $G$, the total system mass $M$, and scale length $\alpha$ are all set equal to 1. The models can then easily be scaled relative to each other by multiplying the initial phase space coordinates and masses by scalar quantities. Furthermore, in these models, each particle represents a large collection of stars since these simulations support $N \propto 10^4$ particles, or many orders of magnitude less than real galaxies. Thus, the greater the number of particles, the higher the simulation resolution. Certain features seen in merger remnants such as stellar shells require large numbers of particles to resolve. As such, simulating millions of particles requires both heavy computational power from super-clusters and more efficient N-body integration schemes.

<figure>
  <img src="ANIMATIONS/panelplot.png" width="1000" align = 'center'>
</figure>

<br>


## Troubleshooting
This code has been optimized for 3-dimensional gravitational interactions of stellar particles. The Numba compiler also expects the input arrays to have the correct dimensions and will error otherwise. For N particles, please ensure the following NumPy arrays have shapes:
- Positions [N x 3]
- Velocities [N x 3]
- Masses [N x 1]
<br>
This is because in 3 dimensions, positions and velocities have x, y, and z components, while the mass array should contain the mass of each particle in the simulation. If one of the arrays has the incorrect shape, please reshape it using the .reshape() NumPy method;
e.i: positions = positions.reshape(N,3) where N is an integer corresponding to the total number of particles.

## Acknowledgments
I would like to thank Professor Jeffrey Kenney, Professor Marla Geha, and Shashank Dattathri for providing invaluable help in completing this project.
This would not have been possible without them. I would also like to thank my Astro Big-Sib Sebastian Monzon and Barry Chiang for helping me out with running the simulations. 

<figure>
  <img src="ANIMATIONS/sim3_grid3x3_dark.png" width="1000" height="1000" align = 'center'>
</figure>


# Documentation
## Table of contents
### [Simulation Setup](#simulation-setup)
### [Running the Simulation](#running-the-simulation)
### [Simulation Analysis](#simulation-analysis)

Here is a demonstration of the functionality of the package in greater detail. A simulation starts with the setup of initial conditions.

```python
from MSG_Nbody import *
```

<a id="simulation-setup"></a>
## Simulation Setup

Initial conditions are loaded into python using the [load_initial_conditions ](https://github.com/elkogerville/MSG_Nbody/blob/main/MSG_Nbody/simulation_setup.py#L14) method. Initial conditions are assumed to be a Nx7 .txt file containing the $x,y,z$ positions, $vx,vy,vz$ velocities and masses $m$ of each particle n. Any initial conditions can be used as long as an Nx3 position array, an Nx3 velocity array, and an Nx1 mass array are provided.

```python
path_2_file = 'initial_conditions/model_disk_3000.txt'
glxy_pos, glxy_vel, glxy_mass = load_initial_conditions(path_2_file)
```

It is often required to manipulate the initial conditions to properly set up a galaxy merger. MSG_Nbody provides a number of functions to facilitate this process. Most importantly, galaxy initial conditions are computed such that the galaxy is in energetic equilibrium when at rest. Thus, great care must be taken when scaling initial conditions. For example, to scale our galaxy's mass and radius, use the [scale_initial_positions](https://github.com/elkogerville/MSG_Nbody/blob/main/MSG_Nbody/simulation_setup.py#L40) method.

```python
# scale galaxy mass and radius by a factor of 10
new_radius = 10
new_mass = 10
glxy_pos, glxy_vel, glxy_mass = scale_initial_positions(glxy_pos, glxy_vel, glxy_mass, new_radius, new_mass)
```

We can also rotate the disk about a specified axis using a rotation matrix, which will rotate the disk around the $x, y,$ or $z$ axis with the [rotate_disk](https://github.com/elkogerville/MSG_Nbody/blob/main/MSG_Nbody/simulation_setup.py#L68) method.

```python
# 45º rotation around y axis
glxy_pos, glxy_vel = rotate_disk(glxy_pos, glxy_vel, 45, 'y')
```

We can compute the escape velocity at a point r=[x,y,z] from the center of mass of the galaxy using the [compute_escape_velocity](https://github.com/elkogerville/MSG_Nbody/blob/main/MSG_Nbody/simulation_setup.py#L112) method.

```python
# escape velocity at point P a distance |P| from a galaxy centered at [0,0,0]
P = [40,40,50]
escape_velocity = compute_escape_velocity(P0[0], P0[1], P0[2], np.sum(glxy_mass))
```

The [concatenate_initial_conditions](https://github.com/elkogerville/MSG_Nbody/blob/main/MSG_Nbody/simulation_setup.py#L139) method allows for the concatenation of an arbitrary number of sets of initial conditions into single ascontiguous positions, velocities, and mass arrays. 
```python
# concatenate 3 sets of galaxy initial conditions, and save the resulting Nx7 initial conditions array to a .txt file
# where N is the sum of the number of particles in each galaxy
pos_list = [glxy1_pos, glxy2_pos, glxy3_pos]
vel_list = [glxy1_vel, glxy2_vel, glxy3_vel]
mass_list = [glxy1_mass, glxy2_mass, glxy3_mass]
positions, velocities, masses = concatenate_initial_conditions(pos_list, vel_list, mass_list, save_2_disk=True)
```

<a id="running-the-simulation"></a> 
## Running the Simulation

To run the simulation, use the [MSG_Nbody](https://github.com/elkogerville/MSG_Nbody/blob/main/MSG_Nbody/MSG_Nbody.py#L20) method. This will create a new folder in the directory the function is ran from and save every 10 timesteps as a Nx7 $x,y,z,vx,vz,vy,\phi$ .npy file, where $\phi$ is the potential each particle feels.
```python
dt = 0.1
timesteps = 2000
MSG_nbody(positions, velocities, masses, dt, timesteps, snapshot_save_rate=10)
```

The gravitational acceleration and potential felt by each particle due to the interactions of each particle is computed using a softened Newtonian potential in the [compute_accel_potential](https://github.com/elkogerville/MSG_Nbody/blob/main/MSG_Nbody/acceleration_potential.py#L16) method.

<a id="simulation-analysis"></a>
## Simulation Analysis

To load simulation snapshots back into python, use the [load_simulation_outputs](https://github.com/elkogerville/MSG_Nbody/blob/main/MSG_Nbody/input_output.py#L21) method. This will separate each timestep into an arbitrary number of subarrays. Each returned object is a list of len(N_particles) containing the separated position, velocity, and potential array of each galaxy. The arrays each have shapes TxNx3 for positions and velocities, and TxNx1 for masses, where T is the number of timesteps.
```python
path_2_snapshots = 'simulation_outputs_N6000/*'
# number of particles per galaxy
N_particles = [3000, 3000]
positions, velocities, potentials = load_simulation_outputs(path_2_snapshots, N_particles)
```

To shift the positions and velocities to a specified frame of reference, use the [shift_2_com_frame](https://github.com/elkogerville/MSG_Nbody/blob/main/MSG_Nbody/analysis.py#L23) method.
```python
# shift all particles to simulation center of mass frame
positions, velocities = shift_2_com_frame(positions, velocities, masses)
# shift all particles to galaxy 1 center of mass frame
# this centers everything around glxy1
positions, velocities = shift_2_com_frame(positions, velocities, gxy1_mass, galaxy_idx=0)
```

To display a simulation snapshot at a timestep $t$, use the [display_galaxies](https://github.com/elkogerville/MSG_Nbody/blob/main/MSG_Nbody/analysis.py#L83) method. Keep in mind that if every 10 timesteps are saved, the total number of snapshots available to plot is the number of timesteps divided by snapshot_save_rate.
```python
timestep = 300
display_galaxies(positions, timestep, sort=True, scale=30)
```

Alternatively, a 3x3 grid plot showing 9 timesteps can be plotted with the [plot_3x3grid](https://github.com/elkogerville/MSG_Nbody/blob/main/MSG_Nbody/analysis.py#L867) method.
```python
# plot x,y projection
t = [50, 230, 260, 300, 400, 500, 700, 800, 900]
plot_grid3x3(positions, t, [0,1], sort=True, snapshot_save_rate=10, savefig=True)
```

We can compute the relative Energy per timestep using the [compute_relative_energy](https://github.com/elkogerville/MSG_Nbody/blob/main/MSG_Nbody/analysis.py#L162) method. This returns a list of containing a TxNx1 energy array for each galaxy.
```python
energies = compute_relative_energy(velocities, potentials)
```

To plot the log distribution of energies for a given galaxy use the [plot_Ne](https://github.com/elkogerville/MSG_Nbody/blob/main/MSG_Nbody/analysis.py#L361) method. A list of timesteps to plot can be passed in.
```python
# plot the log energy distribution of galaxy 1 at timesteps 0, 2600, and 9000
t = [0, 260, 900]
plot_Ne(energies[0], t, snapshot_save_rate=10, savefig=True)
```

To plot a simulated position-velocity diagram (PVD) of an orthagonal projection of a simulation snapshot along a specified line of sight, use the [plot_PVD](https://github.com/elkogerville/MSG_Nbody/blob/main/MSG_Nbody/analysis.py#L463) method.
```python
# PVD of galaxy 2 at timestep 2000 along the z line of sight
timestep = 200
line_of_sight = [0,0,1]
slice_width = 0.4
plot_PVD(positions[1], velocities[1], timestep, line_of_sight, slice_wifth, snapshot_save_rate=10)
```

To plot a simulation timestep with density histogram subplots, use the [plot_density_histogram](https://github.com/elkogerville/MSG_Nbody/blob/main/MSG_Nbody/analysis.py#L867) method.
```python
# xz projection of timestep 0
plot_density_histogram(positions, 0, [0,2], sort=True, scale=55)
```

</div>
