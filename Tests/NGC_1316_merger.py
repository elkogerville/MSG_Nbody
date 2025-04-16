'''
Author: Elko Gerville-Reache
Created: April 17, 2025
Updated: April 17, 2025

10:1 minor merger of a spherical Hernquist galaxy and a disk galaxy.
This simulation investigates a possible merger scenario reproducing 
the elliptical galaxy NGC 1316, thought to have underwent a collision
~1Gyr ago with a smaller spiral galaxy.

'''

from MSG_Nbody import *

# load initial conditions
disk_pos, disk_vel, disk_mass = load_initial_conditions('Initial_Conditions/disk_galaxy_N80000.txt')
sphr_pos, sphr_vel, sphr_mass = load_initial_conditions('Initial_Conditions/sphr_galaxy_N50000.txt')

# scale satellite galaxy 
initial_mass = np.sum(disk_mass)
# divide by the total mass to normalize the satellite galaxy mass
R = 1/initial_mass
M = 1/initial_mass
disk_pos, disk_vel, disk_mass = scale_initial_positions(disk_pos, disk_vel, disk_mass, R, M)
sphr_pos, sphr_vel, sphr_mass = scale_initial_positions(sphr_pos, sphr_vel, sphr_mass, 2, 10)
print(f'final satellite disk galaxy mass: {np.sum(disk_mass):.4}')
print(f'final spherical galaxy mass: {np.sum(sphr_mass):.4}')

# calculate the escape velocity of the satellite galaxy placed at x=50,y=10,z=0
P0 = [50.0, 10.0, 0.0]
escape_velocity = compute_escape_velocity(P0[0], P0[1], P0[2], np.sum(sphr_mass))
print(f'escape velocity: {escape_velocity:.3}')

# set initial velocity to set the disk galaxy on a collision course 
V0 = [-0.5, 0.0, 0.0]
print('magnitude of satellite galaxy trajectory velocity: ', np.sqrt(V0[0]**2 + V0[1]**2 + V0[2]**2))

# move satellite galaxy to its initial position
disk_pos += P0
# give the satellite galaxy a velocity in the negative x direction less than the escape velocity to put it on a collision path
# the velocity is negative because, from the satellite galaxy, the host galaxy is in the negative x direction
disk_vel += V0

# concatenate initial conditions
pos_list = [sphr_pos, disk_pos]
vel_list = [sphr_vel, disk_vel]
mass_list = [sphr_mass, disk_mass]
positions, velocities, masses = concatenate_initial_conditions(pos_list, vel_list, mass_list)

# run N-body simulation
MSG_nbody(positions, velocities, masses, 0.1, 7000, snapshot_save_rate=10)