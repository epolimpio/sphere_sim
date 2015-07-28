import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
#import cPickle as pickle # For Python 3 is just pickle
import pickle


projection=1 # 0 - cartesian projection, 1 - XY projection
infile='sphere_data.p'

data=pickle.load(open( infile, "rb" ) )

r_vs_t=data[0]
F_vs_t=data[1]
n_vs_t=data[2]

# I should also save the parameters for each simulation, that way I can load
# them rather than do this ugly stuff:
N=len(r_vs_t[0])
L=4
R=0.4
rho=np.sqrt(L**2/(4*np.pi))

def draw_circle(r):
    th=np.linspace(0,2*np.pi,100)
    x=r*np.cos(th)
    y=r*np.sin(th)    
    return(x,y)

def get_e_th(r):
    l=np.sqrt(np.sum(r**2))
    th=np.arccos(r[2]/l)
    phi=np.arctan2(r[1],r[0])
    return np.array([np.cos(th)*np.cos(phi),np.cos(th)*np.sin(phi),-np.sin(th)])
    
def get_e_phi(r):
    phi=np.arctan2(r[1],r[0])
    return np.array([-np.sin(phi),np.cos(phi),0])

# setup figure, draw background
def setup_figure():
    fig=plt.figure(1)
    plt.clf()
    if projection==0:
        ax = plt.axes(xlim=(-0.5,2*np.pi+0.5), ylim=(-0.5,np.pi+0.5))
        plt.plot(2*np.pi*np.array([0,1,1,0,0]),np.pi*np.array([0,0,1,1,0]),'--k')
    
        ax.set_aspect('equal')
    elif projection==1:
        ax = plt.axes(xlim=(-(rho+0.5),(rho+0.5)), ylim=(-(rho+0.5),(rho+0.5)))
        (xo,yo)=draw_circle(rho)
        plt.plot(xo,yo,'--k')
        ax.set_aspect('equal')

    cells=[]
    forces=[]
    arrows=[]
    for i in range(0,N):
        cells += ax.plot([], [], '.', color='k')
        forces += ax.plot([], [], '-b',lw=1)
        arrows += ax.plot([], [], '-r', lw=1)
    
    return(fig,cells,forces,arrows)    
#
# initialization function: plot the background of each frame
def init():
    for i in range(0,N):
        cells[i].set_data([], [])
        forces[i].set_data([], [])
        arrows[i].set_data([], [])
    return cells, forces, arrows

# animation function.  This is called sequentially
def animate(f):
    # load data
    r=r_vs_t[f]
    F=F_vs_t[f]
    n=n_vs_t[f]

    # plot cells that are present
    for i in range(0,len(r)):
        if projection==1:
            cells[i].set_data(r[i,0],r[i,1])
            forces[i].set_data([r[i,0],r[i,0]+F[i,0]],[r[i,1],r[i,1]+F[i,1]])
            arrows[i].set_data([r[i,0],r[i,0]+n[i,0]],[r[i,1],r[i,1]+F[i,1]])
        elif projection==0:
            l=np.sqrt(np.sum(r[i,:]**2))
            th=np.arccos(r[i,2]/l)
            phi=np.arctan2(r[i,1],r[i,0])
            if phi<0:
                phi=phi+2*np.pi
            cells[i].set_data(phi,th)

            e_th=get_e_th(r[i,:])    
            e_phi=get_e_phi(r[i,:])
            
            x_phi=np.inner(F[i,:],e_phi)
            x_th=np.inner(F[i,:],e_th)
            forces[i].set_data([phi,phi+x_phi],[th,th+x_th])

            x_phi=np.inner(n[i,:],e_phi)
            x_th=np.inner(n[i,:],e_th)
            arrows[i].set_data([phi,phi+x_phi],[th,th+x_th])
    return cells, forces, arrows
    
(fig,cells,forces,arrows)=setup_figure()

anim = animation.FuncAnimation(fig, animate, init_func=init,frames=len(r_vs_t), interval=1, blit=False)
plt.show()


