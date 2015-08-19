#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
from matplotlib.patches import FancyArrowPatch, Rectangle
#import cPickle as pickle # For Python 3 is just pickle
import pickle
from datetime import datetime
import sys
from myutils import readConfigFile
from metadata_lib import findFilesWithParameters, transformParameters
from os.path import isfile, join
from fortran_lib import simforces, stripack
from sim_plane_fortran import getDelaunayTrianglesOnPlane, startPositions

def openPickleFile(filename):

    try:
        data=pickle.load(open(filename, "rb" ) )
    except:
        print("ERROR!")
        sys.exit("Could not read file. Make sure that the file %s exists." % filename)

    return data

# Read the date time parameters
if len(sys.argv)>=2:
    str_date = sys.argv[1]
    try:
        yy,mm,dd,h,m,s = [int(x) for x in str_date.split(',')[:6]]
    except:
        print("ERROR!")
        sys.exit("Invalid date. The format must be Y,m,d,H,M,S separated by comma.")
    
    video_filename = './data/plane_data_video_%Y-%m-%d_%H-%M-%S.p'
    dtime = datetime(yy,mm,dd,h,m,s)
    outfile=datetime.strftime(dtime,video_filename)
    data = openPickleFile(outfile)
else:
    # get the last file with the parameters in parameters.ini
    parameters = readConfigFile('parameters.ini')
    metadata = './data/metadata_plane.dat'
    files, dates = findFilesWithParameters(metadata, parameters)
    if files:
        dtime = dates[-1]
    else:
        sys.exit("No file found with the parameters in parameters.ini")
    data = openPickleFile(datetime.strftime(dtime,parameters['outfile_video_plane']))

parameters=transformParameters(data[0])
r_vs_t=data[1]
n_vs_t=data[2]
F_vs_t=data[3]
pairs = data[5]
list_, boundary = getDelaunayTrianglesOnPlane(r_vs_t[0][0:2,:])

# parameters
N = parameters['N']
n_steps = parameters['n_steps']
n_save = parameters['n_save']
phi_pack = parameters['phi_pack']

data_analysis = openPickleFile(datetime.strftime(dtime,parameters['outfile_analysis_plane']))
p_angular = np.array(data_analysis[1])[np.arange(n_save-1, n_steps, n_save)]

output_file = datetime.strftime(dtime,'./movies/sim_plane_%Y-%m-%d_%H-%M-%S.mp4')

# Plot Voronoi?
plot_voronoi = False

# Plot Springs?
plot_springs = True

r_ini, Nx, Ny, b = startPositions(N)
N = Nx*Ny
Lx = 2*(Nx-1)
Ly = np.sqrt(3)*(Ny-1)

# setup figure, draw background
def setup_figure():
    fig=plt.figure(1)
    plt.clf()
    ax = fig.add_subplot(1,1,1)
    ax.set_xlim([0,Lx+1])
    ax.set_ylim([0,Ly+1])


    ax.add_patch(Rectangle((0, 0), Lx, Ly, fill=False))

    cells=[]
    springs=[]
    borders=[]
    for i in range(0,N):
        c = plt.Circle((-0,0),0.5,color=cm.copper(0))
        cells.append(ax.add_artist(c))

    if plot_springs:
        for i in range(0,len(pairs)):
            springs += ax.plot([], [], color=cm.spectral(0))

    if plot_voronoi:
        for i in range(0, pairs2.shape[0]):
            borders += ax.plot([], [], color='k')

    return(fig,cells,springs,borders)

# animation function.  This is called sequentially
def animate(f):
    # load data
    F=F_vs_t[f]
    r=r_vs_t[f]
    n=n_vs_t[f]
    
    for i in range(0,N):
        cells[i].center=(r[0,i],r[1,i])

    if plot_springs:
        for i in range(0,len(pairs)):
            i1 = pairs[i,0] - 1
            i2 = pairs[i,1] - 1
            dist = np.sqrt(np.sum((r[:,i1]- r[:,i2])**2))
            c=int((dist-1)*128)
            springs[i].set_data([r[0,i1], r[0,i2]], [r[1,i1], r[1,i2]])
            springs[i].set_color(cm.spectral(c))

    # if plot_voronoi:
    #     list_, baricenters, out_polygon_dict, pairs2, all_areas = getVoronoiOnSphere(r)
    #     b = rho*baricenters
    #     for i in range(0,len(pairs)):
    #         i1 = pairs2[i,0]
    #         i2 = pairs2[i,1]
    #         borders[i].set_data([b[x_plane,i1], b[x_plane,i2]], [b[y_plane,i1], b[y_plane,i2]])      

      
    return (cells,springs,borders)

plt.clf()
(fig,cells,springs,borders)=setup_figure()
anim = animation.FuncAnimation(fig, animate, frames=n_steps//n_save, interval=1, blit=False)
FFMpegWriter = animation.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib',
    comment='Moving cells')
writer = FFMpegWriter(fps=5, metadata=metadata)
anim.save(output_file, writer=writer)
#plt.show()