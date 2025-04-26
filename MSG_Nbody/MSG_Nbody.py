'''
Author: Elko Gerville-Reache
Date Created: 2023-05-20
Date Modified: 2025-04-12
Description:
    function to simulate a group of particles under the influence of gravity
    given initial positions, velocities, and masses
Dependencies:
    - numpy
    - tqdm
'''

import os
import numpy as np
from tqdm import tqdm
from .simulation_setup import ascontiguousarray
from .acceleration_potential import compute_accel_potential
from .input_output import save_snapshot_2_binary_file, create_output_directory

def MSG_Nbody(positions, velocities, masses, dt, timesteps, **kwargs):
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
    kwargs
    ------
    snapshot_save_rate: int
        the frequency (in terms of timesteps) at which simulation
        snapshots are saved. if snapshot_save_rate = 10, a snapshot
        will be saved to disk every 10 timesteps
    start_idx: int
        last computed timestep. allows for resuming an interrupted simulation
        while ensuring that subsequent snapshots are named sequentially to the
        last saved timestep. ex: if start_idx = 2000, the simulation will start
        from timestep 2001. the initial conditions passed into MSG_Nbody should
        be the positions and velocities of the particles at timestep 2000
    directory: str
        path to computation directory. by default is automatically created
        but it can be useful to pass in an existing directory to resume
        a simulation. if resuming a simulation, remember to set a start_idx
    block_size: int
        size of the blocks used for batch processing of interactions.
        if N < block_size, the algorithm degrades to processing the
        interactions in one calculation. by default, block_size = 3000
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
    if (positions.shape != (N, 3)) or \
       (velocities.shape != (N, 3)) or \
       (masses.shape != (N, 1)):
        error_message = (
            "ERROR: Ensure 'positions' and 'velocities' have shapes (N, 3), "
            "and 'masses' has shape (N, 1), where N is the number of particles.\n"
            f"provided shapes -> positions: {positions.shape}, "
            f"velocities: {velocities.shape}, masses: {masses.shape}\n"
            r"/ᐠ_ ꞈ _ᐟ\ <(Fix it...)"
        )
        raise SystemExit(error_message)
    # save every 10 timesteps
    snapshot_save_rate = 10
    if 'snapshot_save_rate' in kwargs:
        snapshot_save_rate = int(kwargs['snapshot_save_rate'])
    # simulation loop start index
    start_idx = 1
    if 'start_idx' in kwargs:
        start_idx = int(kwargs['start_idx'])
    if 'directory' in kwargs and os.path.isdir(kwargs['directory']):
        if start_idx == 1:
            print(
                "\nWARNING: you have not set a value for 'start_idx' "
                'while providing a directory path.\n'
                'WARNING: this can overwrite simulation snapshots!\n'
                'by default, start_idx is set to 1 for a new simulation.\n'
                'if you are trying to resume a simulation\n'
                'set start_idx to the last completed timestep.\n'
                'ex:\n'
                "data = np.load('snapshot_N100000_timestep_4120.npy')\n"
                'pos, vel = data[:,0:3], data[:,3:6]\n'
                'MSG_Nbody(pos, vel, masses, dt, 5000, directory=directory, start_idx=4120)'
            )
            while True:
                confirm = input("\n\ntype 'yes' or 'y' to proceed: ")
                if confirm.lower() in ['y', 'yes']:
                    break
                else:
                    raise ValueError("aborted: confirmation not received")
        directory = kwargs['directory']
    else:
        directory = create_output_directory(N)
    # batch processing block size
    block_size = 3000
    if 'block_size' in kwargs:
        block_size = int(kwargs['block_size'])
    # ensure integer amount of timesteps
    timesteps = int(timesteps)

    # compute softening length based on Dehnen et al., 2001
    softening = 0.017 * ( (N/1e5) )**(-0.23)
    softening_sq = softening**2

    # allocate acceleration matrix and ensure contiguous arrays
    positions, velocities, masses, accel, potential = ascontiguousarray(positions,
                                                                        velocities,
                                                                        masses)
    # calculate initial accelerations
    accel, potential = compute_accel_potential(positions, masses,
                                               accel, potential,
                                               softening_sq, N,
                                               block_size=block_size)
    # save initial conditions if starting from timestep = 0
    sim_snapshot = np.zeros((N, 7))
    if start_idx == 1:
        save_snapshot_2_binary_file(sim_snapshot, positions, velocities,
                                    potential, directory, N, 0)
        np.save(f'masses_N{N}.npy', masses)
    else:
        timesteps += start_idx
        start_idx += 1
    # simulation loop
    # ---------------
    print(r'simulation running....  /ᐠ –ꞈ –ᐟ\<[pls be patient]')
    for i in tqdm(range(start_idx, timesteps+1)):
        # 1/2 kick
        velocities += accel * dt/2.0

        # drift
        positions += velocities * dt

        # update accelerations
        accel, potential = compute_accel_potential(positions, masses,
                                                   accel, potential,
                                                   softening_sq, N,
                                                   block_size=block_size)
        # update velocities
        velocities += accel * dt/2.0

        # write positions, velocities, and potential to binary file
        if i % snapshot_save_rate == 0:
            save_snapshot_2_binary_file(sim_snapshot, positions, velocities,
                                        potential, directory, N, i)

    print(r'simulation complete [yay!!! =＾● ᆺ ●＾= ✿✧･ﾟ: ✧･ﾟ]')
