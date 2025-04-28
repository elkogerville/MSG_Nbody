'''
author: Elko Gerville-Reache
date Created: 2025-02-10
last Modified: 2025-04-24
purpose: animates simulation snapshots where particles are represented as density points using plt.hexbin
usage: change the params in PARAMS, notably the path to the snapshot directory, number of particles, etc..

notes:
- requires the installation of MSG_Nbody, numpy, matplotlib, celluloid, ffmpeg, and tqdm

'''
from MSG_Nbody import *
import matplotlib.pyplot as plt
from celluloid import Camera
from tqdm import tqdm

#######################################################################################

# PLOT PARAMS
user_cmaps = None
dark_mode = False
axes = [0,1]
ax1, ax2 = axes
scale = 100
gridsize = 300
sort = True

# ANIMATION PARAMS
animation_name = 'my_animation.mp4'
animation_length = 10 # [seconds]
dpi = 300

# SIMULATION PARAMS
directory = '/path/to/simulation/directory'
snapshots = 'simulation_outputs_N3000/*'
N_per_galaxy = [3000, 3000]

#######################################################################################
positions, vel, pot = load_simulation_outputs(directory + snapshots, N_per_galaxy)
timesteps = positions[0].shape[0]-1

extent = [-scale,scale,-scale,scale]
fig = plt.figure(figsize=(7,7))
style = 'dark_background' if dark_mode else 'default'
plt.style.use(style)
camera = Camera(fig)
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['font.family'] = 'Courier New'
plt.rcParams['mathtext.default'] = 'regular'
plt.minorticks_on()
plt.tick_params(axis='both', length=2, direction='in',
                which='both', right=True, top=True)
ax1, ax2 = axes
labels = ['X', 'Y', 'Z']
extent = [-scale, scale, -scale, scale]
positions, _, cmaps = set_plot_colors(positions, False,
                                      user_cmaps=user_cmaps,
                                      dark_mode=dark_mode)

plt.xlabel(labels[ax1], size=16)
plt.ylabel(labels[ax2], size=16)
plt.xlim(-scale, scale)
plt.ylim(-scale, scale)
plt.tight_layout()

N = 1 if sort else len(cmaps)
if sort:
    positions = [np.concatenate(positions, axis=1)]
counter = 0
for t in tqdm(range(timesteps)):
    for i, pos in enumerate(positions):
        plt.hexbin(pos[t,:,ax1], pos[t,:,ax2], gridsize=gridsize,
                   bins='log', extent=extent, cmap=cmaps[counter%N])
        counter += 1
    camera.snap()
plt.close()

# generate animation
animation = camera.animate()
animation.save(animation_name, dpi=dpi, fps=timesteps/animation_length)
