import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
#import cPickle as pickle # For Python 3 is just pickle
import pickle

R=1
data=pickle.load(open( './data/sphere_data.p', "rb" ) )

parameters=data[0]
r_vs_t=data[1]
n_vs_t=data[2]
F_vs_t=data[3]

# number of particles
N = int(parameters['N'])
n_steps = int(parameters['n_steps'])
# setup figure, draw background
def setup_figure():
    fig=plt.figure(1)
    plt.clf()

    ax = plt.axes(xlim=(-5,5), ylim=(-5,5))
    ax.set_aspect('equal')

    cells=[]
    forces=[]
    directions=[]
    centers=[]
    for i in range(0,N):
        c = plt.Circle((-0,0),R,color=cm.copper(0))
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
anim = animation.FuncAnimation(fig, animate,frames=n_steps, interval=1, blit=False)
plt.show()