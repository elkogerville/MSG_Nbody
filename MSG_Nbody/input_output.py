'''
Author: Elko Gerville-Reache
Date Created: 2025-03-17
Date Modified: 2025-04-12
Description:
    functions to manage the I/O operations of the simulation. creates a unique
    directory to save snapshots and manages the snapshot file naming
Dependencies:
    - numpy
    - matplotlib
    - tqdm
'''

import os
import glob
import time
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def load_simulation_outputs(directory, N_per_galaxy):
    '''
    Load and organize simulation output files based on the number of galaxies.
    This function searches for all snapshot files in the given directory,
    Sorts them numerically by timestep, and loads particle data for an arbitrary
    amount of galaxies. The data is then separated into position, velocity,
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
        list of TxNx3 arrays of velocities, where T is the number
        of timesteps, N is the number of particles per galaxy,
        and 3 is the number of dimensions
    velocities : list of np.ndarray[np.float64]
        list of numpy arrays, each with shape (T, N, 3) storing the velocity
        components (vx, vy, vz) for each galaxy
    potential : list of np.ndarray[np.float64]
        list of numpy arrays, each with shape (T, N, 1) containing the
        gravitational potential of each particle
    Notes
    -----
    - the function assumes that each snapshot file is a NumPy array with shape
    (total_particles, 7), where columns 0-2 are positions, 3-5 are velocities,
    and 6 is the gravitational potential. by default MSG_Nbody() will save all
    output files with the correct format
    - particles are assigned to galaxies sequentially based on `N_per_galaxy`
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

def save_snapshot_2_binary_file(snapshot_array, pos, vel,
                                potential, directory, N, i):
    '''
    save the current simulation snapshot to a binary file in working directory
    Parameters
    ----------
    snapshot_array: np.ndarray[np.float64]
        Nx7 array to store particle positions, velocities, and potentials
    pos: np.ndarray[np.float64]
        Nx3 array containing the [x, y, z] positions of all particles
    vel: np.ndarray[np.float64]
        Nx3 array containing the [x, y, z] velocities of all particles
    potential: np.ndarray[np.float64]
        Nx1 array containing the potential felt by each particle from
        the combined gravitational influence of all other particles
    directory: str
        path to directory where snapshot will be saved
    i: int
        current simulation timestep
    '''
    snapshot_array[:, 0:3] = pos
    snapshot_array[:, 3:6] = vel
    snapshot_array[:, 6:7] = potential
    filename = os.path.join(directory, f'snapshot_N{N}_timestep_{i}.npy')
    np.save(filename, snapshot_array)

def create_output_directory(N):
    '''
    creates a directory for storing simulation outputs
    if a directory with the specified name already exists, a unique name
    is generated by appending a timestamp. this ensures simulation
    outputs from different runs are stored separately
    Parameters
    ----------
    N : int
        number of particles in the simulation, used to name the directory
    Returns
    -------
    directory_name : str
        name of the created directory, where simulation snapshots will be saved
    Example
    -------
    -> create_output_directory(10000)
    'simulation_outputs_N10000'  # if unique
    'simulation_outputs_N10000_2025-03-18-14-30-55'  # if a conflict occurs
    '''
    directory_name = f'simulation_outputs_N{N}'
    # create a unique directory if one already exists
    if os.path.exists(directory_name):
        timestamp = time.strftime('%Y-%m-%d-%H-%M-%S')
        directory_name = f"{directory_name}_{timestamp}"
    os.makedirs(directory_name)

    return directory_name

def save_figure_2_disk(dpi):
    '''
    Saves current figure to disk, and prompts user for a filename and format
    Parameters
    ----------
    dpi: float
        resolution in dots per inch
    '''
    file_name = input('input filename for image (ex: myimage.pdf): ')
    format = input('please enter format: png or pdf')
    while format not in {'png', 'pdf'}:
        format = input('please enter format: png or pdf')
    plt.savefig(file_name, dpi=dpi, format=format)

def error_handling_axes(axes):
    '''
    Ensures axes parameter is correctly entered in plotting functions
    Parameters
    ----------
    axes: list of int
        list of len(2) specifying which projection to plot where
        0→X, 1→Y, and 2→Z
        example: axes=[0,1] specifies the xy projection
    Returns
    -------
    axes: list of int
        list of len(2) where elements are in range [0,2] and are integers
    '''
    if len(axes) != 2:
        error = ('projection axes must be a list of length 2 \n'
                'ex: [0,1] would specify x,y projection')
        raise ValueError(error)
    if axes[0] not in {0,1,2} or axes[1] not in {0,1,2}:
        raise ValueError("projection axes must contain values from 0,1,2 only")
    if axes[0] == axes[1]:
        raise ValueError("projection axes must contain two distinct axes")
    axes = [int(x) for x in axes]

    return axes

def error_handling_size(s, positions):
    '''
    Ensures scatter size parameter is correctly formatted for plotting functions
    Parameters
    ----------
    s: float
        size of scatter marker
    positions: list of np.ndarray[np.float64]
        list of TxNx3 arrays of positions, where T is the number
        of timesteps, N is the number of particles per galaxy,
        and 3 is the number of dimensions
    Returns
    -------
    s: list of float
        list of scatter marker size for each array in positions
    '''
    if isinstance(s, list):
        s = s if len(s) >= len(positions) else s * len(positions)
    else:
        s = [s]*len(positions)

    return s
