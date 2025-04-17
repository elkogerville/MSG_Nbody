'''
Author: Elko Gerville-Reache
Date Created: 2023-05-20
Date Modified: 2025-04-13
Description:
    function for the computation of the gravitational acceleration and potential
    experienced by a group of particles due to their combined gravitational attraction
                    rⱼ-rᵢ                             mⱼ
     gᵢ = G∑ⱼmⱼ–––––––––––––––––    (1)    ϕᵢ = G∑ⱼ–––––––––    (2)
               [|rⱼ-rᵢ|² + ϵ²]³ᐟ²                  |rⱼ-rᵢ+ϵ²|
Dependencies:
    - numpy
    - numba
'''
import numpy as np
from numba import njit, prange

@njit(parallel=True, fastmath={'nnan', 'ninf'})
def compute_accel_potential(pos, mass, accel, potential,
                            softening_sq, N, block_size=3000):
    '''
    Computes the gravitational acceleration and potential for each particle due
    to all others using softened Newtonian gravity. The interactions are
    processed in blocks of size `block_size` to improve performance
    Parameters
    ----------
    pos: np.ndarray[np.float64]
        Nx3 array containing the [x, y, z] positions of all particles
    mass: np.ndarray[np.float64]
        Nx1 array containing the mass of each particle
    accel: np.ndarray[np.float64]
        Nx3 array to store the computed gravitational acceleration [ax, ay, az]
        for each particle
    potential: np.ndarray[np.float64]
        Nx1 array to store the computed gravitational potential phi for each
        particle
    softening_sq: float
        square of softening length to prevent division by zero and to define the
        simulation resolution
        sqrt(softening_sq) defines the smallest resolvable scale of interaction
    N: int
        number of particles in simulation
    block_size: int, optional
        size of the blocks used for batch processing of interactions.
        if N < block_size, the algorithm degrades to processing the
        interactions in one calculation. by default, block_size = 3000
    Returns
    -------
    accel: np.ndarray[np.float64]
        Nx3 array of the computed gravitational acceleration [ax, ay, az] for
        each particle
    potential: np.ndarray[np.float64]
        Nx1 array of the gravitational potential experienced by each particle
        due to all other particles
    '''
    # gravitational constant
    G = 1
    # zero accel and potential arrays
    accel.fill(0.0)
    potential.fill(0.0)
    # compute number of blocks
    N_blocks = (N + block_size - 1) // block_size
    # loop through each block
    for i in prange(N_blocks):
        # start index of block I
        start_i = block_size*i
        # end index of block I
        end_i = min((i + 1) * block_size, N)
        # size of block
        I = end_i - start_i
        # seperate block I positions into x,y,z components
        x_i = np.ascontiguousarray(pos[start_i:end_i, 0]).reshape(I, 1)
        y_i = np.ascontiguousarray(pos[start_i:end_i, 1]).reshape(I, 1)
        z_i = np.ascontiguousarray(pos[start_i:end_i, 2]).reshape(I, 1)
        # masses of particles in block I
        mass_I = np.ascontiguousarray(mass[start_i:end_i]).reshape(I, 1)
        # loop through upper triangle of blocks
        for j in prange(i, N_blocks):
            # start index of block J
            start_j = block_size*j
            # end index of block J
            end_j = min((j + 1) * block_size, N)
            # size of block J
            J = end_j - start_j
            # seperate block J positions into x,y,z components
            x_j = np.ascontiguousarray(pos[start_j:end_j, 0]).reshape(J, 1)
            y_j = np.ascontiguousarray(pos[start_j:end_j, 1]).reshape(J, 1)
            z_j = np.ascontiguousarray(pos[start_j:end_j, 2]).reshape(J, 1)
            # masses of particles in block J
            mass_J = np.ascontiguousarray(mass[start_j:end_j]).reshape(J, 1)
            # compute JxI particle-particle seperation matrix
            delx = x_j.T - x_i
            dely = y_j.T - y_i
            delz = z_j.T - z_i
            # compute JxI distance matrix
            r = np.sqrt(delx**2 + dely**2 + delz**2 + softening_sq)
            # inverse distance matrix
            inv_r = 1/r
            # if computing acceleration of a block onto itself, zero diagonal
            # these represent the potential of a particle
            # onto itself which is unphysical
            if i == j:
                np.fill_diagonal(inv_r, 0.0)
            # compute inverse distance cubed for acceleration calculation
            inv_r3 = inv_r**3

            # calculate acceleration
            delx_invr3 = delx * inv_r3
            dely_invr3 = dely * inv_r3
            delz_invr3 = delz * inv_r3

            # acceleration on particles in block I due to block J
            accel[start_i:end_i, 0:1] += G * np.dot(delx_invr3, mass_J)
            accel[start_i:end_i, 1:2] += G * np.dot(dely_invr3, mass_J)
            accel[start_i:end_i, 2:3] += G * np.dot(delz_invr3, mass_J)
            # potential on particles in block I due to block J
            potential[start_i:end_i] += G * np.dot(inv_r, mass_J)

            # acceleration on particles in block J due to block I
            if i != j:
                accel[start_j:end_j, 0:1] -= G * np.dot(delx_invr3.T, mass_I)
                accel[start_j:end_j, 1:2] -= G * np.dot(dely_invr3.T, mass_I)
                accel[start_j:end_j, 2:3] -= G * np.dot(delz_invr3.T, mass_I)
                # potential on particles in block J due to block I
                potential[start_j:end_j] -= G * np.dot(inv_r.T, mass_I)

    return accel, potential
