import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
#import cPickle as pickle # For Python 3 is just pickle
import pickle
from datetime import datetime
import sys
from myutils import readConfigFile
from metadata_lib import findFilesWithParameters
from os.path import isfile, join
from fortran_lib import simforces, stripack

def getDelaunayTrianglesOnSphere(r):
    """
    Use STRIPACK to get the Delaunay Triangles List,
    assuming that the points are on a sphere.

    INPUT: r(3,N) -> array with 3 coordinates for N particles
    OUTPUT: tri_list(3,2*N-4) 
    """
    # Get the parameters (we do not assume unit sphere)
    N = r.shape[1]
    mod_r = np.sqrt(np.sum(r**2,0))
    x = r[0,:]/mod_r
    y = r[1,:]/mod_r
    z = r[2,:]/mod_r

    # Calculate the Delaunay triangulation
    dist = np.zeros(N)
    iwk = np.zeros(2*N)
    list_,lptr,lend,lnew,near,next,dist,ier = stripack.trmesh(x,y,z,iwk[0:N],iwk[N:],dist)
    if ier == 0:
        # Construct the list of neighbors
        ltri = np.zeros((3,2*N-4))
        nt, ltri, ier = stripack.trlist2(list_,lptr,lend,ltri)
        if ier == 0:
            # we need to correct it to start at 0 and transpose it
            # to be in conformation with Python
            return ltri.T - 1

    # Case it fails
    return None

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
    dtime = dates[-1]
    data = openPickleFile(datetime.strftime(dtime,parameters['outfile_video']))

parameters=data[0]
r_vs_t=data[1]
n_vs_t=data[2]
F_vs_t=data[3]
print(len(data))
if len(data) > 5:
    spring=True
    pairs = data[5]
else:
    spring = False

# parameters
N = int(parameters['N'])
n_steps = int(parameters['n_steps'])
n_save = int(parameters['n_save'])
phi_pack = float(parameters['phi_pack'])
rho=np.sqrt(N/4/phi_pack) # sphere radius 

# setup figure, draw background
def setup_figure():
    fig=plt.figure(1)
    plt.clf()
    ax = fig.add_subplot(1,1,1)
    ax.set_xlim([-rho-1,rho+1])
    ax.set_ylim([-rho-1,rho+1])
    ax.set_aspect('equal')

    cells=[]
    forces=[]
    directions=[]
    centers=[]
    springs=[]
    for i in range(0,N):
        c = plt.Circle((-0,0),0.5,color=cm.copper(0))
        cells.append(ax.add_artist(c))
        forces += ax.plot([], [], '-m',lw=1)
        directions += ax.plot([], [], '-k',lw=1)
        centers += ax.plot([], [], '.k')
        
    if spring:
        for i in range(0,len(pairs)):
            springs += ax.plot([], [], '-k')
    return(fig,cells,directions,forces,centers, springs)

# animation function.  This is called sequentially
def animate(f):
    # load data
    F=F_vs_t[f]
    r=r_vs_t[f]
    n=n_vs_t[f]

    indsort=np.argsort(r[2,:])
    
    for i in range(0,N):
        j=indsort[i]
        c=int((r[2,j]+1)/2*256)
        cells[i].center=(r[0,j],r[1,j])
        cells[i].set_facecolor(cm.copper(c))

    if spring:
        for i in range(0,len(pairs)):
            i1 = pairs[i,0] - 1
            i2 = pairs[i,1] - 1
            if (r[2,i1] > 0) and (r[2,i2] > 0):
                dist = np.sqrt(np.sum((r[:,i1]- r[:,i1])**2))
                springs[i].set_data([r[0,i1], r[0,i2]], [r[1,i1], r[1,i2]])
            else:
                springs[i].set_data([], [])
      
    return (cells,directions,forces,centers, springs)

plt.clf()
(fig,cells,directions,forces,centers, springs)=setup_figure()
anim = animation.FuncAnimation(fig, animate, frames=n_steps//n_save, interval=1, blit=False)
plt.show()