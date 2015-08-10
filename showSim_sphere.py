import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
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
fpacking = float(parameters['phi_pack'])
rho=np.sqrt(N/4/fpacking) # sphere radius

projection=1 # 0 - cartesian projection, 1 - XY projection

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
        ax = plt.axes(xlim=(-(rho+1),(rho+1)), ylim=(-(rho+1),(rho+1)))
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
    for i in range(0,N):
        if projection==1:
            cells[i].set_data(r[0,i],r[1,i])
            forces[i].set_data([r[0,i],r[0,i]+F[0,i]],[r[1,i],r[1,i]+F[1,i]])
            arrows[i].set_data([r[0,i],r[0,i]+n[0,i]],[r[1,i],r[1,i]+F[1,i]])
        elif projection==0:
            l=np.sqrt(np.sum(r[:,i]**2))
            th=np.arccos(r[2,i]/l)
            phi=np.arctan2(r[1,i],r[0,i])
            if phi<0:
                phi=phi+2*np.pi
            cells[i].set_data(phi,th)

            e_th=get_e_th(r[:,i])    
            e_phi=get_e_phi(r[:,i])
            
            x_phi=np.inner(F[:,i],e_phi)
            x_th=np.inner(F[:,i],e_th)
            forces[i].set_data([phi,phi+x_phi],[th,th+x_th])

            x_phi=np.inner(n[:,i],e_phi)
            x_th=np.inner(n[:,i],e_th)
            arrows[i].set_data([phi,phi+x_phi],[th,th+x_th])
    return cells, forces, arrows
    
(fig,cells,forces,arrows)=setup_figure()

anim = animation.FuncAnimation(fig, animate, init_func=init,frames=n_steps//n_save, interval=1000, blit=False)
plt.show()


