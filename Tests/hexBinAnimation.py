'''
author: Elko Gerville-Reache
date Created: 2025-02-10
last Modified: 2025-02-20
purpose: animates simulation snapshots where particles are represented as density points using plt.hexbin
usage: change the params in ANIMATION PARAMS, notably the path to the snapshot directory, number of particles, etc..

notes:
- requires the installation of numpy, matplotlib, celluloid, ffmpeg, and tqdm

'''
import numpy as np
import matplotlib.pyplot as plt
from celluloid import Camera
from tqdm import tqdm
################################
# FUNCTIONS
################################
def loadData(string, i):
    # load data and return Nx3 position, Nx3 velocity, and Nx1 potential arrays
    data = np.load(string.format(i))

    return data[:, 0:3], data[:, 3:6], data[:, 6:7]

def loadGalaxy(string, timesteps, N, tempResolution):
    # number of snapshots
    Ntimesteps = int( (timesteps/tempResolution) + 1)
    # 3 dimensions
    Ncoord     = 3
    # allocate arrays for data
    position   = np.zeros(shape = (Ntimesteps, N, Ncoord) )
    velocity   = np.zeros(shape = (Ntimesteps, N, Ncoord) )
    potential  = np.zeros(shape = (Ntimesteps, N, 1) )
    # loop through snapshots and load data
    for i in tqdm(range(Ntimesteps)):
        position[i, :, :], velocity[i, :, :], potential[i, :, :] = loadData(string, i*tempResolution)

    return position, velocity, potential

def computeCOM(pos, vel, mass):
    COM  = np.sum(pos*mass, axis = 0)/np.sum(mass)
    COMv = np.sum(vel*mass, axis = 0)/np.sum(mass)

    return COM, COMv

def centerOfMass(pos, vel, mass, timesteps, tempResolution):
    # number of snapshots
    Ntimesteps = int( (timesteps/tempResolution) + 1)
    # 3 dimensions
    Ncoord     = 3
    # allocate arrays for data
    COM  = np.zeros(shape = (Ntimesteps, 1, Ncoord) )
    COMv = np.zeros(shape = (Ntimesteps, 1, Ncoord) )
    # loop through snapshots and compute COM, COMv
    for i in range(0, Ntimesteps, 1):
        COM[i, :, :], COMv[i, :, :] = computeCOM(pos[i, :, :], vel[i, :, :], mass)

    return COM, COMv

def animateHexbin(positions, start, stop, scale, gridRes, cmap = 'cividis'):
    # simulation box bounds
    extentMag = 100
    extentVec = [-extentMag, extentMag, -extentMag, extentMag]
    # figure params
    fig    = plt.figure(figsize = (10, 10))
    camera = Camera(fig) # define celluloid camera
    plt.rcParams['font.family'] = 'Courier New' # set font
    plt.minorticks_on()
    plt.tick_params(axis = 'both', length = 5, direction = 'in', which = 'both', right = True, top = True)
    plt.xlabel('X', size = 13)
    plt.ylabel('Y', size = 13)
    plt.xlim(-scale, scale),plt.ylim(-scale, scale)
    pos = positions
    for i in range(start, stop):
        plt.hexbin(pos[i,:,0], pos[i,:,1], gridsize = gridRes, extent = extentVec, bins = 'log', cmap = cmap)
        camera.snap()
    plt.close()

    return camera

#######################################################################################
# ANIMATION PARAMS: CHANGE THESE VALUES
# number of particles
N = 23000
# number of timesteps
timesteps      = 7000
# number of time steps between each recorded snapshot of the simulation
tempResolution = 10
# path to data directory
path2Files     = 'disk_elliptcal/disk_ellip7_{}.npy'
path2Masses    = 'masses_ellip3.npy'
# animation start frame
animStart      = 0
# animation end frame
animEnd        = int(timesteps/tempResolution + 1)
# animation duration [seconds]
animLength     = 10
# animation name
animName       = 'galaxyHexAnimv2.mp4'
# center of mass frame
comFrame = True

# plot params
# #####################################################################################
# figure scale xlim ∈ [-scale, scale]; ylim ∈ [-scale, scale]
scale = 70
# grid resolution
gridRes = 200

#######################################################################################
# create animation
# load data
pos, vel, pot = loadGalaxy(path2Files, timesteps, N, tempResolution)
mass          = np.load(path2Masses)
# if true, convert to COM frame of reference by subtracting COM from positions
if comFrame == True:
    print('computing center of mass')
    # compute center of mass
    COM, COMv = centerOfMass(pos, vel, mass, timesteps, tempResolution)
    # convert to center of mass frame
    pos -= COM
    vel -= COMv

print('generating animation, sozz pls b patient ₊˚♡₊')
camera = animateHexbin(pos, animStart, animEnd, scale, gridRes)
# compile snapshots
animation = camera.animate() 
animation.save(animName, fps = animEnd/animLength)
# stop the empty plot from displaying
plt.close() 
exit()
