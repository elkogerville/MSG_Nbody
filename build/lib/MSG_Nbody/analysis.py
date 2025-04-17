'''
Author: Elko Gerville-Reache
Date Created: 2025-03-17
Date Modified: 2025-04-12
Description:
    functions to analyze simulation outputs such as plotting snapshots or
    energy distributions
Dependencies:
    - numpy
    - scipy
    - matplotlib
    - tqdm
'''
import numpy as np
from scipy import stats
from scipy.integrate import quad
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
from tqdm import tqdm
from .input_output import save_figure_2_disk, error_handling_axes

def shift_2_com_frame(positions, velocities, mass, galaxy_idx=None):
    '''
    Converts the current frame of reference based on the center of mass
    of a group of particles. Subtracts the center of mass from all particle
    positions and velocities across all timesteps. Is useful to ensure that the
    simulation is centered around the center of mass of the simulation or a
    galaxy.
    Parameters
    ----------
    positions: list of np.ndarray[np.float64]
        list of TxNx3 arrays of positions, where T is the number
        of timesteps, N is the number of particles per galaxy,
        and 3 is the number of dimensions
    velocities: list of np.ndarray[np.float64]
        list of TxNx3 arrays of velocities, where T is the number
        of timesteps, N is the number of particles per galaxy,
        and 3 is the number of dimensions
    mass: np.ndarray[np.float64]
        Mx1 array of particle masses where M is the total number of particles
        in the simulation. if galaxy_idx != None, mass should be a Nx1 array
        corresponding to the masses of a subgroup of particles
    galaxy_idx: int, optional
        set to None by default. in default mode, the function will shift to
        the frame of reference of the center of mass of the entire simulation.
        can also be used to specify a subgroup of particles to compute the
        center of mass with respect to that subgroup. in this case,
        galaxy_idx should equal the index of the subgroup of particles in
        'positions' and mass should correspond to the masses of that subgroup
    Examples
    --------
    pos = [gal_pos1, gal_pos2] # shapes [(N1x3), (N2x3)]
    vel = [gal_vel1, gal_vel2] # shapes [(N1x3), (N2x3)]
    mass # shape Mx1 where M = N1 + N2
    -> shift to global center of mass frame
    pos_shift, vel_shift = shift_2_com_frame(pos, vel, mass)

    -> shift to galaxy 1 center of mass frame (galaxy 1 is index 0 in pos)
    gal_1mass = mass[:N1,1] # only masses in galaxy 1
    pos_shift, vel_shift =  shift_2_com_frame(pos, vel, gal1_mass, galaxy_idx=0)
    '''
    total_mass = np.sum(mass)
    if galaxy_idx == None:
        pos_stack = np.concatenate(positions, axis = 1)
        vel_stack = np.concatenate(velocities, axis = 1)
    else:
        pos_stack = positions[galaxy_idx]
        vel_stack = velocities[galaxy_idx]
    # loop through each timestep
    for i in tqdm(range(positions[0].shape[0]), desc='shifting frame of reference'):
        pos = pos_stack[i]
        vel = vel_stack[i]
        com_pos = np.sum(pos*mass, axis=0)/total_mass
        com_vel = np.sum(vel*mass, axis=0)/total_mass
        # loop through each galaxy
        for j in range(len(positions)):
            positions[j][i] = positions[j][i] - com_pos
            velocities[j][i] = velocities[j][i] - com_vel

    return positions, velocities

def plot_2D(pos, t, axes, scale=50, cmap_dict=None, cb_idx=0,
            user_colors=None, user_cmaps=None, snapshot_save_rate=10,
            savefig=False, dpi=300, dark_mode=False):
    '''
    Plot a 2D projection of a simulation snapshot
    Parameters
    ----------
    pos: list of np.ndarray[np.float64]
        list of TxNx3 arrays of positions, where T is the number
        of timesteps, N is the number of particles per galaxy,
        and 3 is the number of dimensions
    t: int
        timestep to plot. because the simulation only saves a snapshot every
        'snapshot_save_rate', the total number of recorded timesteps is
        timesteps/snapshot_save_rate. thus, a simulation ran for 2000 timesteps
        with snapshot_save_rate = 10 will only have 200 timesteps to plot
    axes: list of int
        list or array of length 2 specifying which two axes
        (0 for x, 1 for y, 2 for z) should be used for the projection.
        ex: axes = [0,1] would specify the xy projection
    scale: float, optional
        defines the half-width of the plotting region. the x and y limits
        will be set to (-scale, scale)
    cmap_dict: dict, optional
        dictionary mapping an integer key (galaxy index number in pos)
        to an array of of shape N, where N is the number of particles in pos.
        used to apply a colormapping to that galaxy
        example: map the x velocities of the first galxaxy in pos at timestep t=0
        t = 0
        galaxy_idx = 0
        dim = 0
        cmap_dict = {galaxy_idx: velocities[galaxy_idx][t,:,dim]}
    cb_idx: int, optional
        the index of which cmap to use for the colorbar. by default is set to 0
        and corresponds to the lowest galaxy_idx in dict (see above). incrementing
        cb_idx by 1 will then select the next galaxy_idx in the cmap dict.
    user_colors: list of str, optional
            allows user to override default colors/cmaps with a user
            specified custom list of matplotlib colors/
    user_cmaps: list of str, optional
        allows user to override default cmaps with a user
        specified custom list of matplotlib cmaps
    snapshot_save_rate: int, optional
        the frequency (in terms of timesteps) at which simulation
        snapshots are saved. is used to convert from timestep index
        to actual simulation timestep. by default is set to 10.
        this should match the value of the simulation snapshot_save_rate
    savefig: boolean, optional
        saves the figure to the working directory if True
    dpi: int, optional
        dpi of saved figure
    dark_mode: boolean, optional
        if True, uses matplotlib dark_background style
    '''
    axes = error_handling_axes(axes)
    t = int(t)
    style = 'dark_background' if dark_mode else 'default'
    ec = (0, 0, 0) if dark_mode else (1, 1, 1)
    fc = (0, 0, 0) if dark_mode else (1, 1, 1)
    with plt.style.context(style):
        with plt.rc_context({
            'axes.linewidth': 0.6,
            'font.family': ['Courier New', 'DejaVu Sans Mono'],
            'mathtext.default': 'regular'
        }):
            figsize = (6,6) if not cmap_dict else (6.3, 5)
            fig, ax = plt.subplots(figsize=figsize)
            ax.minorticks_on()
            ax.tick_params(axis='both', length=2, direction='in',
                           which='both', right=True, top=True)
            # set plot colors
            pos, colors, cmaps = set_plot_colors(pos, False, user_colors=user_colors,
                                                 user_cmaps=user_cmaps, cmap_dict=cmap_dict,
                                                 dark_mode=dark_mode)
            ax1, ax2 = axes
            ax_labels = ['X', 'Y', 'Z']
            if cmap_dict == None:
                cmap_dict = {}
            counter = 0
            for i, galaxy in enumerate(pos):
                colors_i = cmap_dict.get(i, None)
                if colors_i is not None:
                    im = ax.scatter(galaxy[t][:,ax1], galaxy[t][:,ax2], s=0.5,
                                    c=cmap_dict[i], cmap=cmaps[counter%len(cmaps)])
                    if counter == cb_idx:
                        cbar = fig.colorbar(im, ax=ax)
                        cbar.ax.set_ylabel(r'$V_{X}$', size=16)
                    counter += 1
                else:
                    ax.scatter(galaxy[t][:,ax1], galaxy[t][:,ax2],
                               s=0.5, color=colors[i])

            plt.text(scale/1.8, scale/1.2, 't = {}'.format(t*snapshot_save_rate),
                     bbox=dict(boxstyle="round", ec=ec,fc=fc,))

            ax.set_xlim(-scale, scale)
            ax.set_ylim(-scale, scale)
            ax.set_xlabel(ax_labels[ax1], size = 16)
            ax.set_ylabel(ax_labels[ax2], size = 16)
            plt.tight_layout()
            # if savefig is True, save figure to directory
            if savefig:
                save_figure_2_disk(dpi)

            plt.show()

def display_galaxies(positions, timestep, sort=False, scale=100,
                     user_colors=None, savefig=False, dpi=300, dark_mode=False):
    '''
    Plot the xy and xz projections of a simulation snapshot
    Parameters
    ----------
    positions: list of np.ndarray[np.float64]
        list of TxNx3 arrays of positions, where T is the number
        of timesteps, N is the number of particles per galaxy,
        and 3 is the number of dimensions
    timestep: int or Boolean
        timestep to plot. because the simulation only saves a snapshot every
        'snapshot_save_rate', the total number of recorded timesteps is
        timesteps/snapshot_save_rate. thus, a simulation ran for 2000 timesteps
        with snapshot_save_rate = 10 will only have 200 timesteps to plot.
        if False, will allow user to plot initial simulation distribution
        when no timesteps have been calculated
    sort: boolean, optional
        if True, will sort the particles by the axes not used for plotting
        to ensure the that dimension is taken into account when plotting.
        for example: axes = [0,1] and sort=True will sort all particles by their
        z height (dimension 2) and plot particles with the smallest z height
        first, ensuring particles that are 'higher' are shown on top
    scale: float, optional
        defines the half-width of the plotting region. the x and y limits
        will be set to (-scale, scale)
    savefig: boolean, optional
        saves the figure to the working directory if True
    dpi: int, optional
        dpi of saved figure
    dark_mode: boolean, optional
        if True, uses matplotlib dark_background style
    '''
    # if only 1 timestep provided (avoids program from breaking)
    if timestep == None:
        positions = [pos[np.newaxis,:,:] for pos in positions]
        timestep = 0
    else:
        timestep = int(timestep)
    style = 'dark_background' if dark_mode else 'default'
    with plt.style.context(style):
        with plt.rc_context({
            'axes.linewidth': 0.6,
            'font.family': ['Courier New', 'DejaVu Sans Mono'],
            'mathtext.default': 'regular'
        }):
            # setup figure with 2 subplots for xy and xz projections
            fig, ax = plt.subplots(1,2, figsize=(10,5))
            # format axes, minorticks, fonts, and plot colors
            for a in ax:
                a.minorticks_on()
                a.tick_params(axis='both', length=3, direction='in',
                    which='both', right=True, top=True)
                a.set_xlim(-scale, scale)
                a.set_ylim(-scale, scale)
                a.set_xlabel(r'X', size=16)
            ax[0].set_ylabel(r'Y', size=16)
            ax[1].set_ylabel(r'Z', size=16)
            # set plot colors
            positions, colors, _ = set_plot_colors(positions, sort,
                                                   user_colors=user_colors,
                                                   dark_mode=dark_mode)

            # plot each array in the galaxies list
            for i, galaxy in enumerate(positions):
                if sort:
                    for i, proj_axes in enumerate([[0,1], [0,2]]):
                        sorted_pos, c_arr, a_arr, _ = sort_positions(galaxy,
                                                                     timestep,
                                                                     proj_axes,
                                                                     colors)
                        ax[i].scatter(sorted_pos[:,0], sorted_pos[:,i+1],
                                      s=0.4, color=c_arr, alpha=a_arr)
                else:
                    # plot x,y projection
                    ax[0].scatter(galaxy[timestep,:,0], galaxy[timestep,:,1],
                                  s=0.05, color=colors[i], alpha=0.9)
                    # plot x,z projection
                    ax[1].scatter(galaxy[timestep,:,0], galaxy[timestep,:,2],
                                  s=0.05, color=colors[i], alpha=0.9)
            plt.tight_layout()
            # if savefig is True, save figure to directory
            if savefig:
                save_figure_2_disk(dpi)

            plt.show()

def compute_relative_energy(velocities, potentials):
    '''
    Computes the relative Energy, epsilon, based on the velocity and potential
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

def compute_hernquist_Ne(M, a, N_stars):
    '''
    Computes the N[E] profile of an analytical hernquist galaxy. this computes
    the number of particles at a given energy level for each energy level
    Parameters
    ----------
    M: float
        mass of spherical hernquist galaxy
    a: float
        scale length of spherical hernquist galaxy
    N_stars: int
        number of particles in spherical hernquist galaxy
    Returns
    -------
    NE: np.ndarray[np.float64]
        number of particles at each energy bin
    center_bins: np.ndarray[np.float64]
        central bins of energy values in logarithmic space
    '''
    def compute_rho(r, M, a):
        '''
        Computes hernquist density profile ρ(r) for spherically symmetric galaxies
        Parameters
        ----------
        r: np.ndarray[np.float64]
            (N,) shaped array of radii to compute density at
        M: float
            mass of galaxy
        a: float
            scale length of galaxy
        Returns
        -------
        rho: np.ndarray[np.float64]
            density ρ(r) computed at each radii r
        '''
        rho = (M*a) / ( 2*np.pi * r*((r + a)**3) )

        return rho

    def compute_phi(r, M, a):
        '''
        Computes hernquist potential ϕ(r) for a spherically symmetric galaxy
        Parameters
        ----------
        r: np.ndarray[np.float64]
            (N,) shaped array of radii to compute potential at
        M: float
            mass of galaxy
        a: float
            scale length of galaxy
        Returns
        -------
        phi: np.ndarray[np.float64]
            potential ϕ(r) computed at each radii r
        '''
        G = 1
        phi = -( (G * M)/(r + a) )

        return phi

    def interp_d2rdp2(psi_value):
        '''
        Interpolates second derivative of psi
        Parameters
        ----------
        psi_value: float
            value of psi to compute d²ρ/dΨ² at using interpolation
        Returns
        -------
        d2r_dp2_interp: float
            value of d²ρ/dΨ² computed at psi_value
        '''
        # reverse order of arrays to ensure they are monotonically increasing
        d2r_dp2_interp = np.interp(psi_value, psi[::-1], d2r_dp2[::-1])

        return d2r_dp2_interp

    def fE_integrand(psi, e):
        '''
        Computes the eddington inversion integrand at a given Ψ and energy value
        Parameters
        ----------
        psi: float
            value of psi
        e: float
            energy value to compute integrand at
        Returns
        -------
        integrand: float
            value of integrand evaluated at psi and e
        '''
        integrand = ( interp_d2rdp2(psi) ) / ( np.sqrt(e - psi) )

        return integrand

    def f(E):
        '''
        Integrates the distribution function from 0 to a given energy level E
                 1    E   d²ρ    1
        f(E) = ––––– ∫   ––––– ––––– dΨ
                √8π² ⁰    dΨ²  √(E-Ψ)
        scipy.quad will generate psi values in between 0 and E to interpolate
        at to compute d²ρ/dΨ²
        Parameters
        ----------
        E: float
            Energy value acting as upper bound of definite integral
        Returns
        -------
        f_E: float
            value of integral between 0 and E
        '''
        # index 0 in quad to remove errors
        f_E = (1/(np.sqrt(8)*np.pi**2)) * quad(fE_integrand, 0, E, args = (E),
                                               epsrel = 1e-5)[0]
        return f_E

    # initialize log-spaced radii from .001 to 100
    R_grid = np.logspace(-3, 3, num = 5000)
    # calculate density ρ(r)
    rho = compute_rho(R_grid, M, a)
    # compute potential Ψ(r) = -ϕ
    psi = - compute_phi(R_grid, M, a)
    # compute derivative dρ/dΨ
    dr_dp = np.gradient(rho, psi)
    # compute second derivative d²ρ/dΨ²
    d2r_dp2 = np.gradient(dr_dp, psi)
    E_levels = np.logspace(np.log10(psi[-1]), np.log10(psi[0]), 5000)
    # center of each Energy level
    center_bins = (E_levels[1:] + E_levels[:-1])/2
    # number of bins
    N = center_bins.shape[0]

    # initialize storing array
    F_E = np.empty(len(center_bins))
    # loop through each energy level
    for e in tqdm(range(N), desc = 'calculating f(E)'):
        # calculate f[E] for current energy level e
        f_E = f(center_bins[e])
        # store f[E]
        F_E[e] = f_E

    def r_m(E, M, a):
        G = 1
        return (G*M)/E - a
    def integrand(r, e, M, a):
        return np.sqrt(2*(-compute_phi(r, M, a) - e)) * r**2
    def g(e, M, a):
        '''
        Computes the density of states g(E) with scipy.quad integration
        '''
        rm = r_m(e, M, a)
        return 16*np.pi**2 * (quad(integrand, 0.001, rm, args = (e, M, a))[0])

    # initialize storing array
    G_E = np.empty((N))
    # loop through each radius
    for e in tqdm(range(N), desc='calculating g(E)'):
        # calculate g[E] for current energy level e
        g_E = g(center_bins[e], M, a)
        # store g[E]
        G_E[e] = g_E

    de = np.diff(E_levels)
    NE = F_E*G_E*de*N_stars

    return NE, center_bins

def plot_Ne(energy, timesteps, bin_min=-3, bin_max=0.35,
            plot_hernquist=False, grayscale=False, snapshot_save_rate=10,
            savefig=False, dpi=300, dark_mode=False):
    '''
    Plot the energy distribution of particles across different timesteps
    on a log-log plot
    Parameters
    ----------
    energy : np.ndarray[np.float64]
        TxNx1 array containing energy values for each particle where T is the
        number of timesteps, and N is the number of particles
    timestep: list of int
        list of timesteps to plot. because the simulation only saves a snapshot
        every 'snapshot_save_rate', the total number of recorded timesteps is
        timesteps/snapshot_save_rate. thus, a simulation ran for 2000 timesteps
        with snapshot_save_rate = 10 will only have 200 timesteps to plot
    bin_min, bin_max : float, optional
        the minimum and maximum values for the logarithmic binning
        defined as min=10**(bin_min) and max=10**(bin_max)
    plot_hernquist: boolean, optional
        if True, will plot the analytical N(E) curve of a hernquist galaxy,
        using the 'compute_hernquist_Ne' method. will promt the user to input:
        M: float
            mass of hernquist galaxy
        a: float
            scale length of hernquist galaxy
        N: int
            number of particles in hernquist galaxy
    grayscale: boolean, optional
        if True, sets plot color to black
    snapshot_save_rate: int, optional
        the frequency (in terms of timesteps) at which simulation
        snapshots are saved. is used to convert from timestep index
        to actual simulation timestep. by default is set to 10.
        this should match the value of the simulation snapshot_save_rate
    savefig : bool, optional
        if True, prompts user for a filename and saves the figure
    dpi: int, optional
        dpi of saved figure
    dark_mode: boolean, optional
        if True, uses matplotlib dark_background style
    '''
    # ensure timesteps is in list format
    if isinstance(timesteps, int):
        timesteps = [timesteps]
    timesteps = [int(t) for t in timesteps]
    style = 'dark_background' if dark_mode else 'default'
    with plt.style.context(style):
        with plt.rc_context({
            'axes.linewidth': 0.6,
            'font.family': ['Courier New', 'DejaVu Sans Mono'],
            'mathtext.default': 'regular'
        }):
            # Define colors
            if grayscale:
                colors = ['k' for k in range(6)]
            if dark_mode:
                colors = ['w', '#DC267F', '#7b68ee', '#F1A0FB', '#5CCCA1', '#6A5ACD']
            else:
                colors = ['k', '#483D8B', '#DC267F', '#42A27D', '#6A5ACD', '#91B515']
            use_colorbar = len(timesteps) > len(colors)
            # setup figure
            figsize = (7,6) if use_colorbar else (6,6)
            fig, ax = plt.subplots(figsize=figsize)
            plt.minorticks_on()
            plt.tick_params(axis='both', length=5, direction='in',
                            which='both', right=True, top=True)
            ls = ['-', '--', '-.', ':']
            if use_colorbar:
                cmap = cm.rainbow if not grayscale else cm.gray
                norm = mcolors.Normalize(vmin=min(timesteps) * snapshot_save_rate,
                                         vmax=max(timesteps) * snapshot_save_rate)
                # normalize timesteps
                color_mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
                colors = [cmap(norm(t * snapshot_save_rate)) for t in timesteps]
            if plot_hernquist:
                print('Hernquist Galaxy Params:')
                M = float(input('mass of galaxy M: '))
                a = float(input('scale length of galaxy a: '))
                N = float(input('number of particles N: '))
                NE, centers = compute_hernquist_Ne(M, a, N)
                plt.plot(centers, NE/np.max(NE), c='mediumslateblue',
                         label=f'Hernquist model \nM = {M}, a = {a}')
            # plot histogram for each timestep
            for i, (t, color) in enumerate(zip(timesteps, colors)):
                label = None if use_colorbar else f't = {t * snapshot_save_rate}'
                bins = np.logspace(bin_min, bin_max, 65)
                hist, edges = np.histogram(energy[t], bins=bins)
                center = (edges[1:] + edges[:-1]) / 2
                ax.step(center, hist / np.max(hist), color=color,
                        lw=0.6, label=label, ls=ls[i%len(ls)])
            # labels and scales
            ax.set_xlabel('E', size=16)
            ax.set_ylabel('N(E)', size=16)
            ax.set_yscale('log')
            ax.set_xscale('log')
            # add a legend if not using colorbar
            if not use_colorbar:
                ax.legend(loc='upper left')
            else:
                if plot_hernquist == True:
                    ax.legend(loc='upper left')
                cbar = fig.colorbar(color_mapper, ax = ax)
                cbar.set_label("Timesteps", size = 13, rotation=270, labelpad=15)
            plt.tight_layout()
            # if savefig is True, save figure to directory
            if savefig:
                save_figure_2_disk(dpi)

            plt.show()

def plot_PVD(pos, vel, timestep, line_of_sight, width, m_shift=1,
             b_shift=0, transpose=False, snapshot_save_rate=10,
             savefig=False, dpi=300, dark_mode=False):
    '''
    Generate a position-velocity diagram (PVD) for a system of particles by
    projecting particles onto a 2D plane orthogonal to the line-of-sight vector,
    then extracting particles within a window around a best-fit line in the
    projection and plotting their projected position vs. line-of-sight velocity
    Parameters:
    -----------
    pos: np.ndarray[np.float64]
        TxNx3 array containing the x,y,z positions for each particle
        where T is the number of timesteps, and N is the number of
        particles in that galaxy
    vel: np.ndarray[np.float64]
        TxNx3 array containing the vx,vy,vz velocities for each particle
    timestep: int
        timestep to plot. because the simulation only saves a snapshot every
        'snapshot_save_rate', the total number of recorded timesteps is
        timesteps/snapshot_save_rate. thus, a simulation ran for 2000 timesteps
        with snapshot_save_rate = 10 will only have 200 timesteps to plot
    line_of_sight: list of float
        list or array of length 3 representing the line-of-sight direction as a
        vector in 3D space. this is used to project the velocities along the
        line of sight. the vector will automatically be normalized
    width: float
        half-thickness of the extraction window around the best-fit line.
        defines the ± limits for particle selection
    m_shift: float, optional
        multiplies the slope of the best-fit line. by default is set to 1
    b_shift: float, optional
        shifts the best-fit line intercept value by adding b_shift to b.
        by default is set to 0
    transpose: boolean, optional
        if True, chooses the 'y' component of the position projection as the
        major axis to plot against the line of sight velocity
    snapshot_save_rate: int, optional
        the frequency (in terms of timesteps) at which simulation
        snapshots are saved. is used to convert from timestep index
        to actual simulation timestep. by default is set to 10.
        this should match the value of the simulation snapshot_save_rate
    savefig: boolean, optional
        saves the figure to the working directory if True
    dpi: int, optional
        dpi of saved figure
    '''
    def los_to_angles(los_vector):
        '''
        Convert line of sight vector to azimuth and elevation angles in degrees
        Parameters
        ----------
        los_vector array-like of floats
            array or list of length 3 containing the x,y,z components of the
            line of sight vector. the vector will be normalized automatically
        Returns
        -------
        azimuth: float
            azimuthal angle of the line of sight vector in degrees
        elevation: float
            angle between the z axis and the line of sight vector in degrees
        '''
        vx, vy, vz = los_vector
        # compute azimuth and elevation in radians
        azimuth = np.arctan2(vy, vx)
        elevation = np.arcsin(vz / np.linalg.norm(los_vector))

        return np.degrees(azimuth), np.degrees(elevation)

    def project_particles(positions, elev=0, azim=0):
        '''
        Projects 3D particles onto a 2D plane based on elevation and azimuth angles
        to mimic Matplotlib's ax.view_init(elev, azim) functionality.

        Parameters
        ----------
        positions: np.ndarray[np.float64]
            Nx3 array containing the [x, y, z] positions of all particles.
        elev: float, optional
            Elevation angle in degrees (rotation around x-axis).
        azim: float, optional
            Azimuthal angle in degrees (rotation around z-axis).

        Returns
        -------
        projected_pos: np.ndarray[np.float64]
            Nx2 array containing the projected particles.
        '''
        # Convert angles to radians
        elev_rad = np.radians(elev)
        azim_rad = np.radians(azim)

        # Define the line-of-sight vector (Mimicking Matplotlib's camera)
        v = np.array([
            np.cos(elev_rad) * np.cos(azim_rad),
            np.cos(elev_rad) * np.sin(azim_rad),
            np.sin(elev_rad)
        ])

        # Define the "up" vector (default is +z direction)
        up = np.array([0.0, 0.0, 1.0])

        # Create orthonormal basis using Gram-Schmidt
        u = np.cross(up, v)  # Right vector
        if np.linalg.norm(u) < 1e-6:  # If u is too small (e.g., looking straight up/down)
            u = np.array([1.0, 0.0, 0.0])  # Default to x-axis
        u /= np.linalg.norm(u)

        w = np.cross(v, u)  # Up vector
        w /= np.linalg.norm(w)

        # Project particles onto the new basis
        projected_pos = np.column_stack([positions @ u, positions @ w])

        return projected_pos

    def compute_density(x, y):
        '''
        Estimate density of a group of particles with a gaussian
        kde kernel
        Parameters
        ----------
        x: np.ndarray[np.float64]
            Nx1 array representing the x values of each particle
        y: np.ndarray[np.float64]
            Nx1 array representing the y values of each particle
        Returns
        -------
        Z: np.ndarray[np.float64]
            200x200 array of density contours
        '''
        # get min and max
        xmin, xmax = x.min(), x.max()
        ymin, ymax = y.min(), y.max()
        # compute padding
        padding_x = (xmax - xmin) * 0.2
        padding_y = (ymax - ymin) * 0.2
        # generate a 200x200 grid
        xgrid, ygrid = np.mgrid[xmin-padding_x:xmax+padding_x:200j,
                                ymin-padding_y:ymax+padding_y:200j]
        grid = np.vstack([xgrid.ravel(), ygrid.ravel()])
        # stack x and y values into column matrix
        values = np.vstack([x, y])
        # estimate density with gaussian kernel
        kernel = stats.gaussian_kde(values)
        Z = np.reshape(kernel(grid).T, xgrid.shape)

        return xgrid, ygrid, Z

    def compute_kde_limits(xgrid, ygrid, Z, density_threshold=1e-3):
        '''
        Calculate limits based on regions with significant density.
        Parameters
        ----------
        xgrid: np.ndarray[np.float64]
            NxN grid of x values
        ygrid: np.ndarray[np.float64]
            NxN grid of y values
        Z: np.ndarray[np.float64]
            NxN array of density contours
        density_threshold: float, optional
            fraction of maximum density to consider
            (e.g., 0.001 = 0.1% of peak density)
        Returns
        -------
        xmin, xmax: float
            minimum and maxium limits along x axis
        ymin, ymax: float
            minimum and maxium limits along y axis
        '''
        # normalize Z
        Z_norm = Z / Z.max()
        # find grid points above threshold
        inliers = Z_norm > density_threshold
        # compute mins and maxs
        if np.any(inliers):
            xmin = xgrid[inliers].min()
            xmax = xgrid[inliers].max()
            ymin = ygrid[inliers].min()
            ymax = ygrid[inliers].max()
        else:
            # fallback to data range if no points meet threshold
            xmin, xmax = xgrid.min(), xgrid.max()
            ymin, ymax = ygrid.min(), ygrid.max()

        return xmin, xmax, ymin, ymax

    def make_square_limits(xmin, xmax, ymin, ymax):
        '''
        Expands the smaller range to match the larger one,
        ensuring equal aspect ratio
        Parameters
        ----------
        xmin, xmax: float
            minimum and maximum x values in plot
        ymin, ymax: float
            minimum and maximum y values in plot
        Returns
        -------
        xmin_new, xmax_new: float
            shifted minimum and maximum x values in plot
        ymin_new, ymax_new: float
            shifted minimum and maximum y values in plot
        '''
        x_range = xmax - xmin
        y_range = ymax - ymin
        max_range = max(x_range, y_range)

        x_center = (xmin + xmax) / 2
        y_center = (ymin + ymax) / 2

        # Define new equal limits centered at the original center
        xmin_new = x_center - max_range / 2
        xmax_new = x_center + max_range / 2
        ymin_new = y_center - max_range / 2
        ymax_new = y_center + max_range / 2

        return xmin_new, xmax_new, ymin_new, ymax_new

    # index position and velocity by timestep
    pos = pos[timestep]
    vel = vel[timestep]
    N = pos.shape[0]
    # project particles onto plane orthagonal to line of sight
    # vector using Gram-Schmidt orthogonalization
    los_vector = np.asarray(line_of_sight, dtype=np.float64)
    los_vector /= np.linalg.norm(los_vector)
    azim, elev = los_to_angles(los_vector)
    pos_proj = project_particles(pos, elev, azim)
    # compute best fit linear model on position
    linear_fit = np.polyfit(pos_proj[:,0], pos_proj[:,1], 1)
    line_best_fit = np.poly1d(linear_fit)
    # obtain model parameters and shift by ± width
    # this creates 2 parallel lines to mask positions
    # we want the particles that lie in between the line
    m = line_best_fit[1] * m_shift
    b = line_best_fit[0] + b_shift
    # scales width proportionally to slope
    offset = width * np.sqrt(1 + m**2)
    b1, b2 = b + offset, b - offset
    # mask particles
    los_mask = ((pos_proj[:,1] < pos_proj[:,0]*m + b1) &
                (pos_proj[:,1] > pos_proj[:,0]*m + b2))
    pos_los = pos_proj[los_mask]
    # obtain component of velocities along line of sight
    v_los = np.dot(vel[los_mask], los_vector)
    los_label = f"LOS: [{los_vector[0]:.2f}, {los_vector[1]:.2f}, {los_vector[2]:.2f}]"
    # create plot
    style = 'dark_background' if dark_mode else 'default'
    with plt.style.context(style):
        with plt.rc_context({
            'axes.linewidth': 0.6,
            'font.family': ['Courier New', 'DejaVu Sans Mono'],
            'mathtext.default': 'regular'
        }):
            fig, ax = plt.subplots(1, 2, figsize=(13,4),
                                   gridspec_kw={'width_ratios': [1, 3]})
            ax[0].minorticks_on()
            ax[0].tick_params(axis='both', length=2, direction='in',
                              which='both', right=True, top=True)
            ax[1].minorticks_on()
            ax[1].tick_params(axis='both', length=2,
                              direction='in', labelright='on',
                              labelleft=False,  right=True, top=True)
            # plot real space
            # ---------------
            # plot both linear models
            t = timestep*snapshot_save_rate
            # plot entire projection of particles
            ax[0].scatter(pos_proj[:, 0], pos_proj[:, 1], s=0.05, c='#DC267F')
            # plot particles selected in mask
            ax[0].scatter(pos_los[:,0], pos_los[:,1], s=0.05,
                          c='darkslateblue', label=f't = {t}')
            X, Y, Z = compute_density(pos_proj[:, 0], pos_proj[:, 1])
            # get limits
            xmin, xmax, ymin, ymax = compute_kde_limits(X, Y, Z)
            xmin, xmax, ymin, ymax = make_square_limits(xmin, xmax, ymin, ymax)
            # # set plot limits
            ax[0].set_xlim(xmin, xmax)
            ax[0].set_ylim(ymin, ymax)
            x = np.arange(xmin, xmax, 0.1)
            y = m*x + b
            ax[0].plot(x, y+offset, color = 'mediumvioletred')
            ax[0].plot(x, y-offset, color = 'mediumvioletred')
            # plot legend
            ax[0].legend(handlelength=0, handletextpad=0,
                         fancybox=True, loc = 'best')
            ax[0].set_xlabel(fr'$\phi = {azim:.1f}^\circ$', size=16)
            ax[0].set_ylabel(fr'$\theta = {elev:.1f}^\circ$', size=16)
            # plot phase space
            # ----------------
            # plot position - velocity diagram
            ax[1].yaxis.set_label_position('right')
            axis = 0
            if transpose:
                axis = 1
            ax[1].scatter(pos_los[:,axis], v_los, s = .1, c = 'darkslateblue')
            # compute contours with gaussian kde kernel
            X, Y, Z = compute_density(pos_los[:,axis], v_los)
            cs = ax[1].contour(X, Y, Z, levels=np.logspace(-3,10,24),
                               colors='k', alpha=0.6)
            # get limits from contour paths
            if len(cs.allsegs) > 0:
                all_verts = np.vstack([np.vstack(level_segs) for
                                       level_segs in cs.allsegs if level_segs])
                xmin, xmax = all_verts[:,0].min(), all_verts[:,0].max()
                ymin, ymax = all_verts[:,1].min(), all_verts[:,1].max()
                # add 5% padding
                x_pad = 0.05 * (xmax - xmin)
                y_pad = 0.05 * (ymax - ymin)
                xmin, xmax = xmin - x_pad, xmax + x_pad
                ymin, ymax = ymin - y_pad, ymax + y_pad
            else:
                # fallback if contouring failed
                xmin, xmax = X.min(), X.max()
                ymin, ymax = Y.min(), Y.max()
            ax[1].set_xlim(xmin, xmax)
            ax[1].set_ylim(ymin, ymax)
            ax[1].set_xlabel(r'$X_{\perp LOS}$' , size=16)
            ax[1].set_ylabel(r'V$_{LOS}$', size=16)
            plt.title(los_label, size=14)
            plt.tight_layout()
            # if savefig is True, save figure to directory
            if savefig:
                save_figure_2_disk(dpi)

            plt.show()

def plot_density_histogram(positions, timestep, axes, sort=False, scale=100,
                           user_colors=None, savefig=False, dpi=300, dark_mode=False):
    '''
    Plot an orthogonal projection of a timestep with log density histograms
    Parameters
    ----------
    positions: list of np.ndarray[np.float64]
        list of TxNx3 arrays of positions, where T is the number
        of timesteps, N is the number of particles per galaxy,
        and 3 is the number of dimensions
    timestep: int
        timestep to plot. because the simulation only saves a snapshot every
        'snapshot_save_rate', the total number of recorded timesteps is
        timesteps/snapshot_save_rate. thus, a simulation ran for 2000 timesteps
        with snapshot_save_rate = 10 will only have 200 timesteps to plot
    axes: list of int
        list or array of length 2 specifying which two axes
        (0 for x, 1 for y, 2 for z) should be used for the projection.
        ex: axes = [0,1] would specify the xy projection
    sort: boolean, optional
        if True, will sort the particles by the axes not used for plotting
        to ensure the that dimension is taken into account when plotting.
        for example: axes = [0,1] and sort=True will sort all particles by their
        z height (dimension 2) and plot particles with the smallest z height
        first, ensuring particles that are 'higher' are shown on top
    scale: float, optional
        defines the half-width of the plotting region. the x and y limits
        will be set to (-scale, scale)
    savefig: boolean, optional
        saves the figure to the working directory if True
    dpi: int, optional
        dpi of saved figure
    dark_mode: boolean, optional
        if True, uses matplotlib dark_background style
    '''
    axes = error_handling_axes(axes)
    timestep = int(timestep)
    labels = ['X', 'Y', 'Z']
    style = 'dark_background' if dark_mode else 'default'
    with plt.style.context(style):
        with plt.rc_context({
            'axes.linewidth': 0.6,
            'font.family': ['Courier New', 'DejaVu Sans Mono'],
            'mathtext.default': 'regular'
        }):
            fig = plt.figure(figsize=(6, 6))
            # Adjust grid layout to prevent overlap
            gs = fig.add_gridspec(2, 2, width_ratios=(4, 1.2), height_ratios=(1.2, 4),
                                  left=0.15, right=0.9, bottom=0.15, top=0.9,
                                  wspace=0.09, hspace=0.09)  # Reduced spacing
            # create subplots
            ax = fig.add_subplot(gs[1, 0])  # Main scatter plot
            ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)  # Top histogram
            ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)  # Right histogram
            # configure scales and ticks
            ax_histx.set_yscale('log')
            ax_histy.set_xscale('log')
            ax.minorticks_on()
            # tick parameters for main plot
            ax.tick_params(axis='both', length=2, direction='in', which='both',
                           pad=5, right=True, top=True)  # Added pad
            # tick parameters for top histogram (x-axis)
            ax_histx.tick_params(axis='x', direction='in', which='both',
                                 labelbottom=False, bottom=True)
            ax_histx.tick_params(axis='y', direction = 'in', which='both',
                                 left=True, right=True, labelleft=True, pad=5)
            ax_histx.yaxis.set_label_position("left")
            # tick parameters for right histogram (y-axis)
            ax_histy.tick_params(axis='y', direction='in', which='both',
                                 labelleft=False, left=True)
            ax_histy.tick_params(axis='x', direction = 'in', which='both',
                                 bottom=True, top=True, labelbottom=True, pad=5)
            ax_histy.xaxis.set_label_position("bottom")
            # set plot colors
            positions, colors, _ = set_plot_colors(positions, sort,
                                                   user_colors=user_colors,
                                                   dark_mode=dark_mode)

            ax1, ax2 = axes
            # loop through each galaxy
            for i in range(len(positions)):
                if sort:
                    pos_, c, a, c_unique = sort_positions(positions[i], timestep,
                                                          axes, colors)
                    ax.scatter(pos_[:, ax1], pos_[:, ax2], s=0.4, color=c, alpha=a)
                    tags = pos_[:, 3].astype(int)
                    unique_tags = np.sort(np.unique(tags))
                    # loop through each tag and plot histogram
                    for j in range(c_unique.shape[0]):
                        mask = np.where(pos_[:,3] == unique_tags[j])[0]
                        ax_histx.hist(pos_[mask][:, ax1], bins='auto',
                                      color=c_unique[j], histtype='step',
                                      lw=0.8, density=True)
                        ax_histy.hist(pos_[mask][:, ax2], bins='auto',
                                      orientation='horizontal', color=c_unique[j],
                                      histtype='step', lw=0.8, density=True)
                else:
                    pos = positions[i][timestep]
                    # main scatter plot
                    ax.scatter(pos[:, ax1], pos[:, ax2], s=0.4, color=colors[i])
                    # top histogram (x-axis)
                    ax_histx.hist(pos[:, ax1], bins='auto', color=colors[i],
                                  histtype='step', lw=1, density=True)
                    # right histogram (y-axis)
                    ax_histy.hist(pos[:, ax2], bins='auto', orientation='horizontal',
                                  color=colors[i], histtype='step', lw=1, density=True)
            # set axis limits
            ax.set_xlim(-scale, scale)
            ax.set_ylim(-scale, scale)
            # set labels
            ax.set_xlabel(labels[ax1], labelpad=10, size=16)
            ax.set_ylabel(labels[ax2], labelpad=10, size=16)
            # For dynamic axis labels with LaTeX
            ax_histx.set_ylabel(rf'log[N$_{{{labels[ax1]}}}$]', labelpad=10, size=13)
            ax_histy.set_xlabel(rf'log[N$_{{{labels[ax2]}}}$]', labelpad=10, size=13)
            plt.subplots_adjust(left=0.15, right=0.9, bottom=0.15, top=0.9)
            # if savefig is True, save figure to directory
            if savefig:
                save_figure_2_disk(dpi)

            plt.show()

def plot_grid3x3(positions, timesteps, axes, sort=False, scale=50, user_colors=None,
                 snapshot_save_rate=10, savefig=False, dpi=300, dark_mode=False):
    '''
    Plot a 3x3 grid plot of 9 simulation timesteps
    Parameters
    ----------
    positions: list of np.ndarray[np.float64]
        list of TxNx3 arrays of positions, where T is the number
        of timesteps, N is the number of particles per galaxy,
        and 3 is the number of dimensions
    timestep: list of int
        list of timesteps to plot with length 9. because the simulation only
        saves a snapshot every 'snapshot_save_rate', the total number of
        recorded timesteps is timesteps/snapshot_save_rate. thus, a simulation
        ran for 2000 timesteps with snapshot_save_rate = 10 will only have
        200 timesteps to plot
    axes: list of int
        list or array of length 2 specifying which two axes
        (0 for x, 1 for y, 2 for z) should be used for the projection.
        ex: axes = [0,1] would specify the xy projection
    sort: boolean, optional
        if True, will sort the particles by the axes not used for plotting
        to ensure the that dimension is taken into account when plotting.
        for example: axes = [0,1] and sort=True will sort all particles by their
        z height (dimension 2) and plot particles with the smallest z height
        first, ensuring particles that are 'higher' are shown on top
    scale: float, optional
        defines the half-width of the plotting region. the x and y limits
        will be set to (-scale, scale)
    snapshot_save_rate: int, optional
        the frequency (in terms of timesteps) at which simulation
        snapshots are saved. is used to convert from timestep index
        to actual simulation timestep. by default is set to 10.
        this should match the value of the simulation snapshot_save_rate
    savefig: boolean, optional
        saves the figure to the working directory if True
    dpi: int, optional
        dpi of saved figure
    dark_mode: boolean, optional
        if True, uses matplotlib dark_background style
    '''
    # error handling
    axes = error_handling_axes(axes)
    if len(timesteps) != 9:
        error = ('timesteps should be a list of length 9 \n'
                'ex: plot first 9 snapshots: [0,1,2,3,4,5,6,7,8]')
        raise ValueError(error)
    timesteps = [int(t) for t in timesteps]
    labels = ['X', 'Y', 'Z']
    # plot grid
    # ---------
    style = 'dark_background' if dark_mode else 'default'
    ec = (0, 0, 0) if dark_mode else (1, 1, 1)
    fc = (0, 0, 0) if dark_mode else (1, 1, 1)
    with plt.style.context(style):
        with plt.rc_context({
            'axes.linewidth': 0.6,
            'font.family': ['Courier New', 'DejaVu Sans Mono'],
            'mathtext.default': 'regular'
        }):
            fig = plt.figure(figsize = (10,10))
            # define 3x3 grid of subplots
            gs = fig.add_gridspec(3, 3, hspace = 0, wspace = 0)
            (axs) = gs.subplots(sharex=True, sharey=True)
            # store params for each plot
            labeltop = [[True,True,True], [False,False,False],
                        [False,False,False]]
            labelright = [[False,False,True], [False,False,True],
                          [False,False,True]]
            # generate each subplot
            counter = 0
            for i, a in enumerate(axs):
                for j in range(3):
                    a[j].minorticks_on()
                    a[j].tick_params(axis='both', length=1.7, direction='in',
                                     which='both', labeltop=labeltop[i][j],
                                     labelright=labelright[i][j],
                                     right=True, top=True)
                    a[j].xaxis.set_major_locator(ticker.MaxNLocator(3))
                    a[j].yaxis.set_major_locator(ticker.MaxNLocator(3))
                    counter += 1
            # set plot colors
            positions, colors, _ = set_plot_colors(positions, sort,
                                                   user_colors=user_colors,
                                                   dark_mode=dark_mode)
            ax1, ax2 = axes
            axs[0][1].set_xlim(-scale, scale)
            axs[0][1].set_ylim(-scale, scale)
            # # loop through each galaxy
            for i in range(len(positions)):
                counter = 0
                pos = positions[i]
                # loop through each plot
                for j in range(3):
                    # loop through each subplot
                    for k in range(3):
                        timestep_label = timesteps[counter]*snapshot_save_rate
                        timestep = timesteps[counter]
                        # if True, sort particles by height
                        if sort:
                            pos_, c, a, _ = sort_positions(pos, timestep,
                                                           axes, colors)
                            axs[j][k].scatter(pos_[:,ax1], pos_[:,ax2],
                                              s=0.1, color=c, alpha=a)
                        # default grid plot
                        else:
                            pos_t = pos[timestep]
                            axs[j][k].scatter(pos_t[:,ax1], pos_t[:,ax2],
                                              s=0.1, color=colors[i])
                        # timestep legend
                        axs[j][k].text(-scale*0.9, scale*0.9, f't = {timestep_label}',
                                       size=10, bbox=dict(boxstyle="round",
                                                          ec=ec,fc=fc,))
                        counter += 1

            axs[2][1].set_xlabel(labels[ax1], size = 16)
            axs[1][0].set_ylabel(labels[ax2], size = 16)
            plt.tight_layout()
            # if savefig is True, save figure to directory
            if savefig:
                save_figure_2_disk(dpi)

            plt.show()

# def set_plot_colors(positions, sort, cmap='rainbow_r',
#                     trim=False, dark_mode=False):
#     '''
#     Set plot colors based on plotting parameters
#     Parameters
#     ----------
#     positions: list of np.ndarray[np.float64]
#         list of TxNx3 arrays of positions, where T is the number
#         of timesteps, N is the number of particles per galaxy,
#         and 3 is the number of dimensions
#     sort: boolean
#         determines whether or not to concatenate positions into
#         a single array for sorting
#     cmap: matplotlib.pyplot cmap, optional
#         sets the cmap of the plot
#     trim: boolean, optional
#         if True, returns a color array with len(positions)
#     dark_mode: boolean, optional
#         if True, uses a color palette tuned to the matplotlib
#         dark_background style
#     Returns
#     -------
#     positions:
#         original positions array. if sort=True, the positions array
#         is concatenated along axis=1
#         ex: positions = [ ar1, ar2], sort=True
#         ar1.shape: TxNx3
#         ar2.shape: TxMx3
#         returned positions shape: Tx(N+M)x3
#     colors: array like
#         list of colors for the plot
#     '''
#     colors = ['#483D8B', '#DC267F', '#F1A0FB', '#5CCCA1', '#6A5ACD','k']
#     if dark_mode:
#         colors = ['w', '#DC267F', '#7b68ee', '#F1A0FB', '#5CCCA1', '#6A5ACD']
#     # if sort, tag each particle in positions list
#     # then concatenate into single array
#     if sort:
#         positions = [np.concatenate(tag_particles(positions), axis=1)]
#     else:
#         # if more galaxies than colors, generate new list of colors
#         if len(positions) > len(colors):
#             colors = plt.get_cmap(cmap)(np.linspace(0, 1, len(positions)))
#     if trim:
#         colors = np.array(colors[:len(positions)])

#     return positions, colors

def set_plot_colors(positions, sort, user_colors=None, user_cmaps=None,
                    cmap_dict=None, cmap='rainbow_r', dark_mode=False):
    '''
    Set plot colors for each plotting function based on the arguments
    of the plotting function
    Parameters
    ----------
    positions: list of np.ndarray[np.float64]
        list of TxNx3 arrays of positions, where T is the number
        of timesteps, N is the number of particles per galaxy,
        and 3 is the number of dimensions
    sort: boolean
        determines whether or not to concatenate positions into
        a single array for sorting
    user_colors: list of str, optional
        allows user to override default colors/cmaps with a user
        specified custom list of matplotlib colors/
    user_cmaps: list of str, optional
        allows user to override default cmaps with a user
        specified custom list of matplotlib cmaps
    cmap_dict: dictionary, optional
        dictionary mapping an integer key (galaxy index number in pos)
        to an array of of shape N, where N is the number of particles in pos.
        used to apply a colormapping to that galaxy
        example: map the x velocities of the first galxaxy in pos at timestep t=0
        t = 0
        galaxy_idx = 0
        dim = 0
        cmap = {galaxy_idx: velocities[galaxy_idx][t,:,dim]}
    cmap: matplotlib.pyplot cmap, optional
        sets the cmap that colors are drawn from for the plot
        in the case that more subarrays are provided than colors
    dark_mode: boolean, optional
        if True, uses a color palette tuned to the matplotlib
        dark_background style
    Returns
    -------
    positions:
        original positions array. if sort=True, the positions array
        is concatenated along axis=1
        ex: positions = [ ar1, ar2], sort=True
        ar1.shape: TxNx3
        ar2.shape: TxMx3
        returned positions shape: Tx(N+M)x3
    colors: array like
        list of colors for the plot
    cmaps: array like
        list of matplotlib cmaps for the plot
    '''
    # default color list
    colors = ['w', '#DC267F', '#7b68ee', '#F1A0FB', '#5CCCA1', '#6A5ACD'] \
             if dark_mode else ['#483D8B', '#DC267F', '#F1A0FB', '#5CCCA1', '#6A5ACD','k']
    # default cmap list
    cmaps = ['GnBu_r', 'RdPu_r', 'Purples_r', 'cividis',
             'Grays_r', 'Greens_r', 'BuPu_r', 'summer']

    # number of key:value pairs in cmap dictionary
    N_cmap_dict = len(cmap_dict) if cmap_dict is not None else 0
    N_galaxies = len(positions)
    # concatenate positions if sort = True
    if sort:
        positions = [np.concatenate(tag_particles(positions), axis=1)]
        # if sort = True, cmaps cannot be used
        N_cmap_dict = 0
        if user_cmaps is not None:
            user_cmaps = None
            print('WARNING: user_cmap provided with sort=True\n',
                  'cmaps can only be used with sort=False')
    # number of colors needed is total number of galaxies
    # minus number of galaxies mapped to a cmap
    N_colors_needed = N_galaxies - N_cmap_dict

    # ensure enough colors are provided
    if user_colors is not None:
        # ensure enough colors are provided
        if (len(user_colors) >= N_colors_needed):
            colors = user_colors
        else:
            # if enough default colors to plot all galaxies
            if (N_galaxies <= len(colors)):
                print('WARNING: not enough user specified colors\n',
                      f'number of user colors: {len(user_colors)}\n',
                      f'number of colors needed: {N_colors_needed}\n',
                      f'defaulting to MSG_Nbody colors list: {colors[:N_colors_needed]}')
            # if more galaxies than default colors
            else:
                print('WARNING: not enough user specified colors\n',
                      f'number of user colors: {len(user_colors)}\n',
                      f'number of colors needed: {N_colors_needed}\n',
                      f'defaulting to drawing color sequence from cmap {cmap}')
    # if more galaxies than default colors
    if N_galaxies > len(colors):
        colors = plt.get_cmap(cmap)(np.linspace(0, 1, N_galaxies))

    if user_cmaps is not None:
        # ensure enough cmaps are provided
        if (len(user_cmaps) >= N_cmap_dict):
            cmaps = user_cmaps
        else:
            print('WARNING: not enough user specified cmaps\n',
                  f'number of user cmaps: {len(user_cmaps)}\n',
                  f'number of cmaps needed: {N_cmap_dict}\n',
                  f'defaulting to MSG_Nbody cmaps list: {cmaps[:len(positions)]}')

    return positions, colors, cmaps

def tag_particles(positions):
    '''
    Give each galaxy a unique tag
    Parameters
    ----------
    positions: list of np.ndarray[np.float64]
        list of TxNx3 arrays of positions, where T is the number
        of timesteps, N is the number of particles per galaxy,
        and 3 is the number of dimensions
    Returns
    -------
    pos: list of np.ndarray[np.float64]
        list of TxNx4 arrays of positions, where the new 4th column
        is a unique tag for each galaxy. Thus all particles in each galaxy N
        have a tag of N-1
    '''
    pos = []
    for i in range(len(positions)):
        # number of timesteps
        timesteps = positions[i].shape[0]
        # number of particles
        N = positions[i].shape[1]
        # TxNx1 matrix of ones
        ones = np.ones((timesteps, N, 1))
        # append matrix as a unique tag
        pos_sort = np.append(positions[i], ones*i, axis=2)
        pos.append(pos_sort)

    return pos

def sort_positions(positions, timestep, axes, color_arr):
    '''
    Sort positions along axis omitted in projection and create an array
    of colors with a color corresponding to each particle. For example,
    if plotting xy projection, will sort by z height to ensure that the
    particles are plotted such that higher particles are on top of lower
    ones. an array of colors with 1 color per particle is generated,
    respecting the color of each galaxy
    Parameters
    ----------
    positions: list of np.ndarray[np.float64]
        list of len(1) containing a TxNx4 array of x,y,z positions, where the
        4th column contains a unique tag for each galaxy. Thus all particles
        in each galaxy N have a tag of N-1
    timestep:
        timestep to plot
    axes: list of int
        list or array of length 2 specifying which two axes
        (0 for x, 1 for y, 2 for z) should be used for the projection.
        ex: axes = [0,1] would specify the xy projection
    color_arr: list of str
        list of MSG_Nbody base matplotlib colors
    Returns
    -------
    pos_sorted: np.ndarray[np.float64]
        Nx4 x,y,z,tag positions array sorted by height
    colors: list or array like
        list of colors or RBGA array of colors to use in plot
    alphas: np.ndarray[float]
        array of alpha values for each particle
    unique_colors: np.ndarray[str or float]
        array of each unique color
    '''
    # compute sorting axis (axis not used)
    sorting_axis = 3 - sum(axes)
    # get particles at timestep
    pos = positions[timestep]
    # sort by sorting axis
    sorted_indeces = np.argsort(pos[:, sorting_axis])
    pos_sorted = pos[sorted_indeces]
    # array of all tags sorted by height
    tags = pos_sorted[:, 3].astype(int)
    # get unique tags (1 unique tag per galaxy)
    unique_tags = np.sort(np.unique(tags))
    # generate increasing alphas per galaxy
    min_alpha, max_alpha = 0.6, 0.8
    alpha_step = (max_alpha - min_alpha) / max(1, len(unique_tags)-1)
    # create alpha mapping where each galaxy gets slightly higher alpha
    alpha_map = {tag: min_alpha + i*alpha_step
                for i, tag in enumerate(unique_tags)}
    # apply alphas to all particles
    alphas = np.array([alpha_map[tag] for tag in tags])
    # array to store unique colors
    unique_colors = []
    # if more galaxies than base colors
    # ---------------------------------
    if len(unique_tags) > len(color_arr):
        # generate 1 color per galaxy from cmap
        cmap = plt.get_cmap('rainbow_r', len(unique_tags))
        # generate RGBA color array
        colors = np.zeros((len(tags), 4))
        # for each unique tag assign color to all galaxy particles
        for i, tag in enumerate(unique_tags):
            mask = tags == tag
            colors[mask] = cmap(i)
            unique_colors.append(cmap(i))
        # reshape unique colors
        unique_colors = np.asarray(unique_colors).reshape(-1,4)
    # if less galaxies than base colors
    # ---------------------------------
    else:
        # color dictionary mapping tags to a color
        color_map = {tag: color_arr[i] for i, tag in enumerate(unique_tags)}
        # for each particle assign color mapped to its galaxy
        colors = np.array([color_map[tag] for tag in tags])
        # get unique colors
        for i in color_map:
            unique_colors.append(color_map[i])
        unique_colors = np.asarray(unique_colors)

    return pos_sorted, colors, alphas, unique_colors
