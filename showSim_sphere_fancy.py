import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
#import cPickle as pickle # For Python 3 is just pickle
import pickle
from datetime import datetime
import sys
from myutils import readConfigFile
from metadata_lib import findFilesWithParameters, transformParameters
from os.path import isfile, join
from fortran_lib import simforces, stripack
from sim_sphere_fortran import getDelaunayTrianglesOnSphere

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
    
    video_filename = './data/sphere_data_video_%Y-%m-%d_%H-%M-%S.p'
    dtime = datetime(yy,mm,dd,h,m,s)
    outfile=datetime.strftime(dtime,analysis_filename)
    data = openPickleFile(outfile)
else:
    # get the last file with the parameters in parameters.ini
    parameters = readConfigFile('parameters.ini')
    metadata = './data/metadata.dat'
    files, dates = findFilesWithParameters(metadata, parameters)
    if files:
        dtime = dates[-1]
    else:
        sys.exit("No file found with the parameters in parameters.ini")
    data = openPickleFile(datetime.strftime(dtime,parameters['outfile_video']))

parameters=transformParameters(data[0])
r_vs_t=data[1]
n_vs_t=data[2]
F_vs_t=data[3]
pairs = data[5]

# parameters
N = parameters['N']
n_steps = parameters['n_steps']
n_save = parameters['n_save']
phi_pack = parameters['phi_pack']
update_nn = parameters['update_nn'] == 1
rho=np.sqrt(N/4/phi_pack) # sphere radius 

# setup the plane of the plot
x_plane = 0
y_plane = 2
z_plane = 1
# if True, then invert z
z_inverted = False
if z_inverted:
    for i in range(0,len(r_vs_t)):
        r_vs_t[i][z_plane,:] = -r_vs_t[i][z_plane,:]

# setup figure, draw background
def setup_figure():
    fig=plt.figure(1)
    plt.clf()
    ax = fig.add_subplot(1,1,1)
    ax.set_xlim([-rho-1,rho+1])
    ax.set_ylim([-rho-1,rho+1])
    ax.set_aspect('equal')

    cells=[]
    springs=[]
    for i in range(0,N):
        c = plt.Circle((-0,0),0.5,color=cm.copper(0))
        cells.append(ax.add_artist(c))

    for i in range(0,len(pairs)):
        springs += ax.plot([], [], color=cm.spectral(0))

    return(fig,cells,springs)

# animation function.  This is called sequentially
def animate(f):
    global pairs
    # load data
    F=F_vs_t[f]
    r=r_vs_t[f]
    n=n_vs_t[f]
    if update_nn:
        pairs = simforces.get_all_pairs(getDelaunayTrianglesOnSphere(r)+1)
    

    indsort=np.argsort(r[z_plane,:])
    
    for i in range(0,N):
        j=indsort[i]
        c=int((r[z_plane,j]+1)/2*256)
        cells[i].center=(r[x_plane,j],r[y_plane,j])
        cells[i].set_facecolor(cm.copper(c))

    for i in range(0,len(pairs)):
        i1 = pairs[i,0] - 1
        i2 = pairs[i,1] - 1
        if (r[z_plane,i1] > 0) and (r[z_plane,i2] > 0):
            dist = np.sqrt(np.sum((r[:,i1]- r[:,i2])**2))
            c=int((dist-1)*128)
            springs[i].set_data([r[x_plane,i1], r[x_plane,i2]], [r[y_plane,i1], r[y_plane,i2]])
            springs[i].set_color(cm.spectral(c))
        else:
            springs[i].set_data([], [])
      
    return (cells,springs)

plt.clf()
(fig,cells,springs)=setup_figure()
anim = animation.FuncAnimation(fig, animate, frames=n_steps//n_save, interval=1, blit=False)
plt.show()