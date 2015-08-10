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

    ax = plt.axes(xlim=(-rho-1,rho+1), ylim=(-rho-1,rho+1))
    ax.set_aspect('equal')

    cells=[]
    forces=[]
    directions=[]
    centers=[]
    for i in range(0,N):
        c = plt.Circle((-0,0),0.1,color=cm.copper(0))
        cells.append(ax.add_artist(c))
        forces += ax.plot([], [], '-m',lw=1)
        directions += ax.plot([], [], '-k',lw=1)
        centers += ax.plot([], [], '.k')
    
    return(fig,cells,directions,forces,centers)

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
        
    return (cells,directions,forces,centers)

plt.clf()
(fig,cells,directions,forces,centers)=setup_figure()
anim = animation.FuncAnimation(fig, animate, frames=n_steps//n_save, interval=1, blit=False)
plt.show()