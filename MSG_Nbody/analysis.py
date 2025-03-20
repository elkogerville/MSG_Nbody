'''
Author: Elko Gerville-Reache
Date Created: 2025-03-17
Date Modified: 2025-03-19
Description:
    functions to analyze simulation outputs such as plotting snapshots or
    energy distributions
Dependencies:
    - numpy
    - matplotlib
'''
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from tqdm import tqdm

def load_simulation_outputs(directory, N_per_galaxy):
    '''
    load and organize simulation output files based on the number of galaxies.
    this function searches for all snapshot files in the given directory,
    sorts them numerically by timestep, and loads particle data for an arbitrary
    amount of galaxies. the data is then separated into position, velocity,
    and potential arrays for each galaxy
    Parameters
    ----------
    directory : str
        path to the snapshot files (e.g., 'simulation_outputs/*')
        the path to the directory must be given, followed by a '*' which tells
        the glob python module to extract all files in the directory.
        the files must be named in the format:
        'name_XXX.npy', where XXX is an integer timestep
    N_per_galaxy : list of int
        list where each element is the number of particles in a given galaxy
    Returns
    -------
    positions : list of np.ndarray[np.float64]
        A list of numpy arrays, each with shape (T, N, 3), where:
        - T is the number of timesteps
        - N is the number of particles in that galaxy
        - 3 corresponds to the (x, y, z) coordinates
    velocities : list of np.ndarray[np.float64]
        A list of numpy arrays, each with shape (T, N, 3) storing the velocity
        components (vx, vy, vz) for each galaxy
    potential : list of np.ndarray[np.float64]
        A list of numpy arrays, each with shape (T, N, 1) containing the
        gravitational potential of each particle
    Notes
    -----
    - The function assumes that each snapshot file is a NumPy array with shape
    (total_particles, 7), where columns 0-2 are positions, 3-5 are velocities,
    and 6 is the gravitational potential. by default MSG_Nbody() will save all
    output files with the correct format
    - Particles are assigned to galaxies sequentially based on `N_per_galaxy`
    Example
    -------
    -> load a 2 galaxy simulation where galaxy1 has N=12000, and galaxy2 N=9000
     # search for all files inside directory 'simulation_outputs_6000'
    directory = '/simulation_outputs_20000/*'
    N_per_galaxy = [12000, 9000]
    '''
    # search for all files in directory and sort by timestep
    files = sorted(glob.glob(directory),
                   key=lambda x: int(x.split('_')[-1].split('.')[0]))
    # number of timesteps
    timesteps = len(files)
    # number of dimensions (x,y,z)
    N_coordinates = 3
    # number of galaxies to separate
    N_galaxies = len(N_per_galaxy)
    # allocate a set of arrays per galaxy
    positions = [np.zeros((timesteps, N, N_coordinates)) for N in N_per_galaxy]
    velocities = [np.zeros((timesteps, N, N_coordinates)) for N in N_per_galaxy]
    potential = [np.zeros((timesteps, N, 1)) for N in N_per_galaxy]
    # loop through each timestep and load data
    for i in tqdm(range(timesteps)):
        snapshot = files[i]
        # load timestep data
        data = np.load(snapshot)
        pos = data[:, 0:3]
        vel = data[:, 3:6]
        pot = data[:, 6:7]
        # split data into galaxies
        start = 0
        for j, Ngal in enumerate(N_per_galaxy):
            positions[j][i] = pos[start : start + Ngal]
            velocities[j][i] = vel[start : start + Ngal]
            potential[j][i] = pot[start : start + Ngal]
            # move start index for next galaxy
            start += Ngal

    return positions, velocities, potential

def display_galaxies(galaxies, timestep, scale=100, savefig=False):
    '''
    plot the xy and xz projections of a simulation snapshot
    Parameters
    ----------
    galaxies: list of np.ndarray[np.float64]
        list containing each set of particle positions
    timestep: int
        timestep to plot. because the simulation only saves a snapshot
        every 'snapshot_save_rate', the total number of recorded
        timesteps is timesteps/snapshot_save_rate. thus, if a simulation
        is ran for 2000 timesteps and snapshot_save_rate = 10, there
        will only be 200 timesteps to plot
    scale: float
        defines the half-width of the plotting region. the x and y limits
        will be set to (-scale, scale)
    savefig: boolean
        saves the figure to the working directory if True
    '''
    timestep = int(timestep)
    # setup figure with 2 subplots for xy and xz projections
    fig, ax = plt.subplots(1,2, figsize=(10,5))
    # format axes, minorticks, fonts, and plot colors
    plt.rcParams.update({
        'axes.linewidth': 0.6,
        'font.family': 'Courier New',
        'mathtext.default': 'regular'
    })
    for a in ax:
        a.minorticks_on()
        a.tick_params(axis='both', length=3, direction='in',
            which='both', right=True, top=True)
        a.set_xlim(-scale, scale)
        a.set_ylim(-scale, scale)
        a.set_xlabel(r'X', size=15)
    ax[0].set_ylabel(r'Y', size=15)
    ax[1].set_ylabel(r'Z', size=15)
    colors = ['darkslateblue', 'mediumvioletred', 'violet',
              'mediumslateblue', 'orchid', 'k']
    # if more galaxies than defined colors, generate a new list of colors
    if len(galaxies) > len(colors):
        cmap = plt.get_cmap('rainbow')
        colors = [cmap(i) for i in np.linspace(0, 1, len(galaxies))]

    # plot each array in the galaxies list
    for i, g in enumerate(galaxies):
        # plot x,y projection
        ax[0].scatter(g[timestep,:,0], g[timestep,:,1],
                      s = 0.1, color = colors[i])
        # plot x,z projection
        ax[1].scatter(g[timestep,:,0], g[timestep,:,2],
                      s = 0.1, color = colors[i])

    plt.tight_layout()
    # if savefig is True, save figure to directory
    if savefig is not False:
        # promt user to input file name
        file_name = input('please input filename for image (ex: myimage.png): ')
        plt.savefig(file_name, dpi = 300, format = 'png')

    plt.show()

def compute_relative_energy(velocities, potentials):
    '''
    computes the relative Energy, epsilon, based on the velocity and potential
    Parameters
    ----------
    velocity: list of np.ndarray[np.float64]
        list of TxNx3 arrays of velocities, where T is the number
        of timesteps, N is the number of particles per galaxy,
        and 3 is the number of dimensions
    potential: list of np.ndarray[np.float64]
        list of TxNx1 arrays of potential values, where T is the
        number of timesteps, and N is the number of particles
        per galaxy
    Returns
    -------
    epsilons: list of np.ndarray[np.float64]
        list of TxNx1 arrays of relative Energy levels where
        T is the number of timesteps, and N is the number of
        particles per galaxy
    '''
    epsilons = []
    # loop through each set of galaxies
    for velocity, potential in zip(velocities, potentials):
        # compute kinetic energy
        xvel = velocity[:, :, 0:1]
        yvel = velocity[:, :, 1:2]
        zvel = velocity[:, :, 2:3]
        kinetic_energy = (1/2)*(xvel**2 + yvel**2 + zvel**2)
        # compute relative energy
        epsilon = potential - kinetic_energy
        epsilons.append(epsilon)

    return epsilons

def plot_Ne(energy, timesteps, savefig=False, bin_m=-3,
            bin_M=0.35, snapshot_save_rate=10):
    '''
    plot the energy distribution of particles across different timesteps
    Parameters
    ----------
    energy : np.ndarray[np.float64]
        TxNx1 array containing energy values for each particle where T is the number
        of timesteps, and N is the number of particles
    timesteps : list of int
        list of timesteps to plot
    savefig : bool, optional
        if True, prompts user for a filename and saves the figure
    bin_m, bin_M : float, optional
        the minimum and maximum values for the logarithmic binning
    snapshot_save_rate : int, optional
        conversion factor from timestep index to actual simulation time
    '''
    # Define colors
    colors = ['k', 'darkslateblue', 'mediumvioletred', 'violet',
              'mediumslateblue', 'orchid', 'purple']
    use_colorbar = len(timesteps) > len(colors)
    # setup figure
    if use_colorbar:
        fig, ax = plt.subplots(figsize = (7, 6))
    else:
        fig, ax = plt.subplots(figsize = (6, 6))
    plt.rcParams.update({
        'axes.linewidth': 0.6,
        'font.family': 'Courier New',
        'mathtext.default': 'regular'
    })
    plt.minorticks_on()
    plt.tick_params(axis='both', length=5, direction='in',
                    which='both', right=True, top=True)

    if use_colorbar:
        cmap = cm.rainbow
        norm = mcolors.Normalize(vmin=min(timesteps) * snapshot_save_rate,
                                 vmax=max(timesteps) * snapshot_save_rate)
        # normalize timesteps
        color_mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
        colors = [cmap(norm(t * snapshot_save_rate)) for t in timesteps]

    # plot histogram for each timestep
    for i, (t, color) in enumerate(zip(timesteps, colors)):
        label = None if use_colorbar else f't = {t * snapshot_save_rate}'
        bins = np.logspace(bin_m, bin_M, 65)
        hist, edges = np.histogram(energy[t], bins=bins)
        center = (edges[1:] + edges[:-1]) / 2
        ax.step(center, hist / np.max(hist), color=color, lw=0.6, label=label)

    # labels and scales
    ax.set_xlabel('E', size=13)
    ax.set_ylabel('N(E)', size=13)
    ax.set_yscale('log')
    ax.set_xscale('log')

    # add a legend if not using colorbar
    if not use_colorbar:
        ax.legend(loc='upper left')
    else:
        cbar = fig.colorbar(color_mapper, ax = ax)
        cbar.set_label("Timestep", size = 13, rotation=270, labelpad=15)
    plt.tight_layout()
    # save figure if needed
    if savefig:
        file_name = input('Please input filename for image (ex: myimage.png): ')
        plt.savefig(file_name, dpi=300, format='png')

    plt.show()
