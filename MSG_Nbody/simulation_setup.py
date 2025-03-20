'''
Author: Elko Gerville-Reache
Date Created: 2025-03-17
Date Modified: 2025-03-19
Description:
    functions to handle simulation setup, such as loading initial conditions,
    rotating particle positions and velocities, and merging arrays
Dependencies:
    - numpy
'''

import numpy as np

def load_initial_conditions(filename):
    '''
    load initial conditions from a .txt file and separate
    into position, velocity, and mass arrays
    Parameters
    ----------
    filename: str
        path to data, should be a Nx7 txt file
    Returns
    -------
    positions: np.ndarray[np.float64]
        Nx3 array containing the [x, y, z] positions of all particles
    velocities: np.ndarray[np.float64]
        Nx3 array containing the [x, y, z] velocities of all particles
    masses: np.ndarray[np.float64]
        Nx1 array containing the mass of each particle
    '''
    # read initial conditions into Python
    data = np.loadtxt(filename)
    # index the initial N x 7 array into respective arrays
    positions = data[:, 0:3].astype(np.float64)
    velocities = data[:, 3:6].astype(np.float64)
    masses = data[:, 6:7].astype(np.float64)

    return positions, velocities, masses

def scale_initial_positions(pos, vel, mass, R, M):
    '''
    correctly scales initial conditions by a given radius R, and mass M
    Parameters
    ----------
    pos: np.ndarray[np.float64]
        Nx3 array containing the [x, y, z] positions of all particles
    vel: np.ndarray[np.float64]
        Nx3 array containing the [x, y, z] velocities of all particles
    mass: np.ndarray[np.float64]
        Nx1 array containing the mass of each particle
    R: float
        radius to scale positions to
    M: float
        scales masses by this parameter
    Returns
    -------
    returns original arrays scaled by R and M
    '''
    G = 1
    # scale position, velocity, and mass by scalar quantities
    # the velocities are scaled proportionally to a circular orbit velocity
    pos = pos * R
    vel = vel * np.sqrt(G*M/R)
    mass = mass * M

    return pos, vel, mass

def rotate_disk(pos, vel, deg, axis):
    '''
    rotates positions and velocities of a galaxy model about a specified axis
    Parameters
    ----------
    pos: np.ndarray[np.float64]
        Nx3 array containing the [x, y, z] positions of all particles
    vel: np.ndarray[np.float64]
        Nx3 array containing the [x, y, z] velocities of all particles
    deg: float
        degrees to rotate by
    axis: str
        set this argument to 'x', 'y', or 'z' to rotate around that axis
    Returns
    -------
    rotated_pos: np.ndarray[np.float64]
        Nx3 array containing the rotated [x, y, z] positions of all particles
    rotated_vel: np.ndarray[np.float64]
        Nx3 array containing the rotated [x, y, z] velocities of all particles
    '''
    if axis.lower() not in ['x', 'y', 'z']:
        raise ValueError("invalid axis. Choose 'x', 'y', or 'z'!")
    # ensure position and velocity arrays are appropriate shapes for dot product
    N = pos.shape[0]
    if (pos.shape != (N,3)) or (vel.shape != (N,3)):
        error_message = 'ERROR: ensure input arrays have shape Nx3'
        raise ValueError(error_message)
    # convert degrees to radians
    rad = np.deg2rad(deg)
    sin = np.sin(rad)
    cos = np.cos(rad)
    # define x, y, z rotation matrices
    rotation_matrices = {
        'x': np.array([[1, 0, 0], [0, cos, -sin], [0, sin, cos]]),
        'y': np.array([[cos, 0, sin], [0, 1, 0], [-sin, 0, cos]]),
        'z': np.array([[cos, -sin, 0], [sin, cos, 0], [0, 0, 1]])
    }
    # rotate positions and velocities with a dot product
    rotation_matrix = rotation_matrices[axis.lower()]
    rotated_pos = pos @ rotation_matrix
    rotated_vel = vel @ rotation_matrix

    return rotated_pos, rotated_vel

def compute_escape_velocity(x, y, z, M):
    '''
    computes the escape velocity at a given [x, y, z] distance
    from a point mass centered at the origin
    Parameters
    ----------
    x: float
        x distance from point mass
    y: float
        y distance from point mass
    z: float
        z distance from point mass
    M: float
        mass of point mass
    Returns
    -------
    escape_vel: [float]
        escape velocity at a given distance from a point mass with mass M
    '''
    G = 1
    # find length of vector between origin and position
    r = np.sqrt(x**2 + y**2 + z**2)
    # compute escape velocity
    escape_vel = np.sqrt(2*G*M/r)

    return escape_vel

def concatenate_initial_conditions(pos_list, vel_list, mass_list, save_2_disk = False):
    '''
    concatenates the intial conditions of an arbitrary amount of galaxies
    Parameters
    ----------
    pos_list: list of np.ndarray[np.float64]
        list containing each set of particle positions
    vel_list: list of np.ndarray[np.float64]
        list containing each set of particle velocities
    mass_list: list of np.ndarray[np.float64]
        list containing each set of particle masses
    Returns
    -------
    positions : np.ndarray[np.float64]
        an array of shape (sum(N), 3), where sum(N) is the total number of particles
        across all input sets. contains the [x, y, z] positions of all particles
    velocities : np.ndarray[np.float64]
        An array of shape (sum(N), 3). contains the [x, y, z] velocities of all particles
    masses : np.ndarray[np.float64]
        An array of shape (sum(N), 1). contains the masses of all particles
    Example
    -------
    3 galaxies with N1 = 3000, N2 = 4000, & N3 = 5000
        pos_list = [galaxy1_pos, galaxy2_pos, galaxy3_pos]
        vel_list = [galaxy1_vel, galaxy2_vel, galaxy3_vel]
        mass_list = [galaxy1_mass, galaxy2_mass, galaxy3_mass]

    The resulting arrays will have shapes:
        positions   -> (12000, 3)
        velocities  -> (12000, 3)
        masses      -> (12000, 1)
    '''
    # ensure an equal amount of positions, velocities, and masses
    if not (len(pos_list) == len(vel_list) == len(mass_list)):
        raise ValueError('ensure an equal amount of initial conditions per list')
    # ensure proper shapes across each set of initial conditions
    for pos, vel, mass in zip(pos_list, vel_list, mass_list):
        N = pos.shape[0]
        if (pos.shape != (N, 3)) or (vel.shape != (N, 3)) or (mass.shape != (N, 1)):
            error_message = (
            "ERROR: ensure 'pos' and 'vel' have shape (N, 3), "
            "and 'mass' has shape (N, 1), where N is the number of particles \n"
            f"provided shapes -> pos: {pos.shape}, gal_vel: {vel.shape}, mass: {mass.shape}\n"
            r"/ᐠ_ ꞈ _ᐟ\ <(Fix it...)"
        )
            raise ValueError(error_message)

    # concatenate initial conditions
    pos = np.concatenate(pos_list)
    vel = np.concatenate(vel_list)
    mass = np.concatenate(mass_list)
    # create contiguous arrays
    positions = np.ascontiguousarray(pos)
    velocities = np.ascontiguousarray(vel)
    masses = np.ascontiguousarray(mass)
    print(
        f'positions shape: {positions.shape}, velocities shape: {velocities.shape}, '
        f'masses shape: {masses.shape} total simulation mass: {np.sum(masses)}'
        )
    if save_2_disk == True:
        N = positions.shape[0]
        initial_conditions = np.hstack([positions, velocities, masses])
        np.savetxt(f'init_conditions_N{N}.txt', initial_conditions, newline='\n')

    return positions, velocities, masses

def ascontiguousarray(positions, velocities, masses):
    '''
    converts input arrays to contiguous arrays in memory and allocates an
    acceleration array initialized to zeros
    Parameters
    ----------
    positions: np.ndarray[np.float64]
        Nx3 array containing the [x, y, z] positions of all particles
    velocities: np.ndarray[np.float64]
        Nx3 array containing the [x, y, z] velocities of all particles
    masses: np.ndarray[np.float64]
        Nx1 array containing the masses of all particles
    Returns
    -------
    positions: np.ndarray[np.float64]
        contiguous Nx3 array of particle positions
    velocities: np.ndarray[np.float64]
        contiguous Nx3 array of particle velocities
    masses: np.ndarray[np.float64]
        contiguous Nx1 array of particle masses
    accel: np.ndarray[np.float64]
        contiguous Nx3 array of zeros to store particle accelerations
    '''
    N = positions.shape[0]
    positions = np.ascontiguousarray(positions)
    velocities = np.ascontiguousarray(velocities)
    masses = np.ascontiguousarray(masses)
    accel = np.ascontiguousarray(np.zeros((N, 3)))

    return positions, velocities, masses, accel
