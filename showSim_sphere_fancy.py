import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
#import cPickle as pickle # For Python 3 is just pickle
import pickle

R=0.25
data=pickle.load(open( 'sphere_data.p', "rb" ) )

r_vs_t=data[0]
n_vs_t=data[1]
F_vs_t=data[2]

N=len(r_vs_t[0])
# setup figure, draw background
def setup_figure():
    fig=plt.figure(1)
    plt.clf()

    ax = plt.axes(xlim=(-2,2), ylim=(-2,2))
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

    indsort=np.argsort(r[:,2])
    
    for i in range(0,N):
        j=indsort[i]
        c=int((r[j,2]+1)/2*256)
        cells[i].center=(r[j,0],r[j,1])
        cells[i].set_facecolor(cm.copper(c))
        
    return (cells,directions,forces,centers)

plt.clf()
(fig,cells,directions,forces,centers)=setup_figure()
anim = animation.FuncAnimation(fig, animate,frames=len(r_vs_t), interval=1, blit=False)
plt.show()