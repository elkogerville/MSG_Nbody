'''
MSG_Nbody test program

Elko Gerville-Reache
Created: April 12, 2025
Updated: April 12, 2025

-------------------------------------------------------------------------
Run a 10:1 merger simulation between a spherical galaxy and a disk galaxy
with 16,000 total particles

** NOTE **
Please download the initial conditions from the Initial_Conditions folder
in the github repo. Please place these files in the same directory containing
MSG_Nbody.py
'''

from MSG_Nbody import *

# load initial conditions
sphr_file = 'Initial_Conditions/sphr_galaxy_N10000.txt'
disk_file = 'Initial_Conditions/disk_galaxy_N6000'
sphr_pos, sphr_vel, sphr_mass = load_initial_conditions(sphr_file)
disk_pos, disk_vel, disk_mass = load_initial_conditions(disk_file)

# scale disk galaxy
initial_mass = np.sum(disk_mass)
# divide by the total mass to normalize the satellite galaxy mass
R = 1/initial_mass
M = 1/initial_mass
disk_pos, disk_vel, disk_mass = scale_initial_positions(disk_pos, disk_vel,
                                                        disk_mass, R, M)

# scale spherical galaxy initial conditions
sphr_pos, sphr_vel, sphr_mass = scale_initial_positions(sphr_pos, sphr_vel,
                                                        sphr_mass, 2, 10)

print(f'final satellite disk galaxy mass: {np.sum(disk_mass)}')
print(f'final host spherical galaxy mass: {np.sum(sphr_mass)}')

# calculate the escape velocity of the satellite galaxy placed at x=50,y=30,z=15
P0 = [50.0, 10.0, 0.0]
escape_velocity = compute_escape_velocity(P0[0], P0[1], P0[2], np.sum(sphr_mass))
print(f'escape velocity: {escape_velocity}')

V0 = [-0.5, 0.0, 0.0]
ve_magnitude = np.sqrt(V0[0]**2 + V0[1]**2 + V0[2]**2)
print('magnitude of satellite galaxy trajectory velocity: ', ve_magnitude)
# move satellite galaxy to its initial position
disk_pos += P0
# give the satellite galaxy a velocity in the negative x direction less than
# the escape velocity to put it on a collision path
# the velocity is negative because, from the satellite galaxy,
# the host galaxy is in the negative x direction
disk_vel += V0

pos_list = [sphr_pos, disk_pos]
vel_list = [sphr_vel, disk_vel]
mass_list = [sphr_mass, disk_mass]
positions, velocities, masses = concatenate_initial_conditions(pos_list,
                                                               vel_list,
                                                               mass_list)

# display the initial simulation setup!
dt = 0.1
timesteps = 5000
plot_orbital_trajectory(pos_list, vel_list, mass_list, dt,
                        timesteps, scale=80, plot_glxys=True)
plt.close()

# run N-body simulation
MSG_Nbody(positions, velocities, masses, 0.1, 2000, snapshot_save_rate=10)
