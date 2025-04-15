'''
Author: Elko Gerville-Reache
Created: April 14, 2025
Updated: April 14, 2025

Program to check the progress of an ongoing simulation
Please provide the path to the directory containing the simulation snapshots
as well as the number of particles per galaxy under USER INPUTS
'''


from MSG_Nbody import *

# ––––––––––––––––––––––––––––––– USER INPUTS ––––––––––––––––––––––––––––––––––

# enter path to simulation output directory
# e.i 'simulation_outputs_N20000/*'
path_2_simulation_outputs = 'simulation_outputs_N20000/*'
# enter number of particles per galaxy
# e.i [2000, 2000] for 2 galaxies with 2000 particles each
N_per_galaxy = [10000, 10000]
pos, vel, pot = load_simulation_outputs(path_2_simulation_outputs, N_per_galaxy)

# ––––––––––––––––––––––––––––––– PLOT PARAMS ––––––––––––––––––––––––––––––––––
# plot xy projection by default
axes = [0,1]
# plotting range of [-50, 50] in both dimensions
scale = 50

# –––––––––––––––––––––––––––––––––– PLOT ––––––––––––––––––––––––––––––––––––––

# plot simulation
N_timesteps = pos[0].shape[0]
t = np.linspace(0, N_timesteps-1, 9)
plot_grid3x3(pos, t, axes=axes, sort=True, scale=scale)
plt.close()

# plot energies
energies = compute_relative_energy(vel, pot)
t = np.linspace(0, N_timesteps-1, 3)
for i in range(len(energies)):
    plot_Ne(energies[i], t, bin_min=-3, bin_max=3)
    plt.close()
