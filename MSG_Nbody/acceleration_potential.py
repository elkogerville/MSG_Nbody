'''
Author: Elko Gerville-Reache
Date Created: 2023-05-20
Date Modified: 2025-04-12
Description:
    function for the computation of the gravitational acceleration and potential
    experienced by a group of particles
                       rⱼ-rᵢ
     gᵢ = G ∑ⱼ mⱼ –––––––––––––––– 
                   [|rⱼ-rᵢ|² + ϵ²]
Dependencies:
    - numpy
    - numba
'''
import numpy as np
from numba import njit

@njit(parallel=True, fastmath={'nnan', 'ninf'})
def compute_accel_potential(pos, mass, accel, potential, softening_sq, N):
    '''
    Computes the gravitational acceleration and potential for each particle due
    to all others using softened Newtonian gravity
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
    Returns
    -------
    accel: np.ndarray[np.float64]
        Nx3 array of the computed gravitational acceleration [ax, ay, az] for
        each particle
    potential: np.ndarray[np.float64]
        Nx1 array of the gravitational potential experienced by each particle
        due to all other particles
    '''
    G = 1
    # seperate position into x,y,z components
    x = np.ascontiguousarray(pos[:, 0]).reshape(N, 1)
    y = np.ascontiguousarray(pos[:, 1]).reshape(N, 1)
    z = np.ascontiguousarray(pos[:, 2]).reshape(N, 1)
    # calculate particle-particle seperations
    delx = x.T - x
    dely = y.T - y
    delz = z.T - z
    r = np.sqrt(delx**2 + dely**2 + delz**2 + softening_sq)
    inv_r3 = r**(-3)
    # calculate acceleration
    accel[:, 0:1] = G * np.dot((delx * inv_r3), mass)
    accel[:, 1:2] = G * np.dot((dely * inv_r3), mass)
    accel[:, 2:3] = G * np.dot((delz * inv_r3), mass)
    accel = np.ascontiguousarray(accel)

    # calculate (N x N) particle-particle potential matrix
    potential_N = (G * mass.reshape(1, N)) / r
    # set diagonal elements to zero
    # these represent the potential of a particle onto itself which is unphysical
    np.fill_diagonal(potential_N, 0.0)
    potential[:] = np.sum(potential_N, axis=0).reshape(N, 1)

    return accel, potential
