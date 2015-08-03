import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
#import cPickle as pickle # For Python 3 is just pickle
import pickle
from datetime import datetime
import sys

if len(sys.argv)>=2:
    str_date = sys.argv[1]
else:
    str_date = input("Insert date, format Y,m,d,H,M,S ->")

try:
    yy,mm,dd,h,m,s = [int(x) for x in str_date.split(',')[:6]]
except:
    print("ERROR!")
    sys.exit("Invalid date. The format must be Y,m,d,H,M,S separated by comma.")

outfile_video=datetime.strftime(datetime(yy,mm,dd,h,m,s),'./data/sphere_data_video_%Y-%m-%d_%H-%M-%S.p')

try:
    data=pickle.load(open(outfile_video, "rb" ) )
except:
    print("ERROR!")
    sys.exit("Could not read file. Make sure that the file %s exists." % outfile_video)

parameters=data[0]
r_vs_t=data[1]
n_vs_t=data[2]
F_vs_t=data[3]

# parameters
N = int(parameters['N'])
n_steps = int(parameters['n_steps'])
n_save = int(parameters['n_save'])
fpacking = float(parameters['fpacking'])
rho=np.sqrt(N/4/fpacking) # sphere radius 

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
        c = plt.Circle((-0,0),1,color=cm.copper(0))
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