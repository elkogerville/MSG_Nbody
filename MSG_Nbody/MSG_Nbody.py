'''
Author: Elko Gerville-Reache
Date Created: 2023-05-20
Date Modified: 2025-03-19
Description:
    function to simulate a group of particles under the influence of gravity
    given initial positions, velocities, and masses
Dependencies:
    - os
    - numpy
    - tqdm
'''

import os
import numpy as np
from tqdm import tqdm
from .simulation_setup import ascontiguousarray
from .acceleration_potential import compute_accel_potential
from .output import save_snapshot_2_binary_file, create_output_directory

def MSG_nbody(positions, velocities, masses, dt, timesteps, **kwargs):
    '''
    runs an N-body simulation to model the dynamics of a system of particles
    under mutual gravitational attraction. the function integrates the equations
    of motion of the particles and advances them using a leap frog integration
    scheme. it saves snapshots of the system's state (positions, velocities,
    and potential) to disk at user-specified intervals. the I/O operations
    are managed automatically by the os package and require no user input
    Parameters
    ----------
    positions: np.ndarray[np.float64]
        Nx3 array containing the [x, y, z] positions of all particles
    velocities: np.ndarray[np.float64]
        Nx3 array containing the [x, y, z] velocities of all particles
    masses: np.ndarray[np.float64]
        Nx1 array containing the masses of all particles
    dt: float
        time step size for advancing the simulation. smaller values reduce
        integration errors but increase runtime
    timesteps: int
        number of timesteps to simulate
    kwargs: dict
        optional keyword arguments:
        - 'snapshot_save_rate': int, optional (default: 10)
            the frequency (in terms of timesteps) at which to save the
            simulation snapshots. The simulation will save a snapshot every
            'snapshot_save_rate' timesteps
    Returns
    -------
    None
        this function outputs a binary file to disk every 'snapshot_save_rate'
        iterations
    '''
    # simulation setup
    # ----------------
    # compute number of particles
    N = positions.shape[0]
    # enforce appropriate array dimensions
    if (positions.shape != (N, 3)) or (velocities.shape != (N, 3)) or (masses.shape != (N, 1)):
        error_message = (
            "ERROR: Ensure 'gal_pos' and 'gal_vel' have shape (N, 3), "
            "and 'mass' has shape (N, 1), where N is the number of particles.\n"
            f"provided shapes -> positions: {positions.shape}, velocities: {velocities.shape}, "
            f"masses: {masses.shape}\n"
            r"/ᐠ_ ꞈ _ᐟ\ <(Fix it...)"
        )
        raise SystemExit(error_message)
    # save every 10 timesteps
    snapshot_save_rate = 10
    if 'snapshot_save_rate' in kwargs:
        snapshot_save_rate = int(kwargs['snapshot_save_rate'])
    # ensure integer amount of timesteps
    timesteps = int(timesteps)
    directory = create_output_directory(N)

    # compute softening length based on Dehnen et al., 2001
    softening = 0.017 * ( (N/1e5) )**(-0.23)
    softening_sq = softening**2

    # allocate acceleration matrix and ensure contiguous arrays
    positions, velocities, masses, accel = ascontiguousarray(positions, velocities, masses)

    # calculate initial accelerations
    accel, potential = compute_accel_potential(positions, masses, accel, softening_sq, N)
    # save initial conditions
    sim_snapshot = np.zeros((N, 7))
    save_snapshot_2_binary_file(sim_snapshot, positions, velocities, potential, directory, N, 0)

    # simulation loop
    # ---------------
    print(r'simulation running....  /ᐠ –ꞈ –ᐟ\<[pls be patient]')
    for i in tqdm(range(1, timesteps+1)):
        # 1/2 kick
        velocities += accel * dt/2.0

        # drift
        positions += velocities * dt

        # update accelerations
        accel, potential = compute_accel_potential(positions, masses, accel, softening_sq, N)

        # update velocities
        velocities += accel * dt/2.0

        # write positions, velocities, and potential to binary file
        if i % snapshot_save_rate == 0:
            save_snapshot_2_binary_file(sim_snapshot, positions, velocities, potential, directory, N, i)

    return r'simulation complete [yay!!! (ﾐΦ ﻌ Φﾐ)✿ *ᵖᵘʳʳ*]'
