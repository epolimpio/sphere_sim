import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri
from scipy.spatial import Delaunay
#import cPickle as pickle # For Python 3 is just pickle
import pickle
from myutils import readConfigFile 
from fortran_lib import simforces
import time
start_time = time.time()

# function for getting unit vectors in spherical coordinates for each position <r>
def get_e_r(r):
    l=np.sqrt(np.sum(r**2))
    return r/l

def get_e_th(r):
    l=np.sqrt(np.sum(r**2))
    th=np.arccos(r[2]/l)
    phi=np.arctan2(r[1],r[0])
    return np.array([np.cos(th)*np.cos(phi),np.cos(th)*np.sin(phi),-np.sin(th)])
    
def get_e_phi(r):
    phi=np.arctan2(r[1],r[0])
    return np.array([-np.sin(phi),np.cos(phi),0])

# function for changing position <r> so that it is on a sphere of radius <rho>
def constrain_to_sphere(r):
    l=np.sqrt(np.sum(r**2))
    return r/l*rho

# ----- PARAMETERS ----- #
parameters = readConfigFile('parameters.ini')

# number of particles
N = int(parameters['N'])

# simulation parameters
F_th_0 = float(parameters['F_th_0'])
mu_th = float(parameters['mu_th'])
mu_F = float(parameters['mu_F'])
eta_n = float(parameters['eta_n'])
K = float(parameters['K'])
R = float(parameters['R'])
n_steps = int(parameters['n_steps'])

# timestep
dt = float(parameters['dt'])

# outfile for data
outfile = str(parameters['outfile'])

# calculate radius of sphere... kind of ugly due to random parameter <L>
L = float(parameters['L'])
rho=np.sqrt(N*L**2/(4*np.pi)) # sphere radius

# ----- INITIALIZE RANDOM POSITION AND ORIENTATION ----- #
# get random positions (i.e. <th>,<phi> in spherical coordinates) for all particles
th=np.pi*np.random.rand(N)
phi=2*np.pi*np.random.rand(N)
# get random angle for unit vector <n> in the local plane
nu=2*np.pi*np.random.rand(N)

r=np.zeros((3,N))
n=np.zeros((3,N))
e_r=np.zeros((3,N))

# calculate positions <r_i> and unit vectors <e_r,i> and <n_i>
for i in range(0,N):
    r[:,i]=rho*np.array([np.sin(th[i])*np.cos(phi[i]),np.sin(th[i])*np.sin(phi[i]),np.cos(th[i])])
    e_r[:,i]=get_e_r(r[:,i])
    n[:,i]=np.cos(nu[i])*get_e_th(r[:,i]) + np.sin(nu[i])*get_e_phi(r[:,i])

# clear variables for storing data    
r_vs_t=[np.array([]),]*n_steps
F_vs_t=[np.array([]),]*n_steps
n_vs_t=[np.array([]),]*n_steps

# ----- RUN FOR SEVERAL TIMES ----- #
for t in range(0,n_steps):
    print(t)

    # caculate total of active and repulsive forces
    F_tot=simforces.calc_force_elastic(r,n,F_th_0,K,R)
    F_tot_plane=np.zeros((3,N))
    # calculate force in the local plane at position <r_i>
    for i in range(0,N):    
        F_tot_plane[:,i]=F_tot[:,i]-np.inner(F_tot[:,i],e_r[:,i])*e_r[:,i]
    
    # integrate equation of motion for direction n_i
    for i in range(0,N):
        # first calculate unit vector of total force
        f=F_tot_plane[:,i]/np.sqrt(np.sum(F_tot_plane[:,i]**2))
        # next, calculate axis of rotation
        e_rot=np.cross(n[:,i],f)
        e_rot_L=np.sqrt(np.sum(e_rot**2))
        if e_rot_L>1e-4:
            e_rot=e_rot/e_rot_L
        # then, calculate absolute angle A between n and f
        B=np.inner(n[:,i],f)
        if B>1.0:
            A=0.0
        elif B<-1.0:
            A=np.pi
        else:             
            A=np.arccos(np.inner(n[:,i],f))
        # finally, calculate new direction n
        dn=dt*(mu_th*A + eta_n*np.random.randn())*np.cross(e_rot,n[:,i])
        n[:,i] = n[:,i] + dn
        n[:,i]=n[:,i]/np.sqrt(np.sum(n[:,i]**2))

    # integrate equation of motion for position <r_i>
    for i in range(0,N):  
        dr = dt*mu_F*F_tot_plane[:,i]
        r[:,i] = r[:,i] + dr        
        r[:,i] = constrain_to_sphere(r[:,i])

    # calculate unit vector <n_i> in the plane at the new position <r_i(t+dt)>
    for i in range(0,N):
        e_r[:,i]=get_e_r(r[:,i])
        n[:,i] = n[:,i] - np.inner(n[:,i],e_r[:,i])*e_r[:,i]
        n[:,i]=n[:,i]/np.sqrt(np.sum(n[:,i]**2))

    # store data for <r_i>, <N-i> and <F_tot,i>
    r_vs_t[t] = np.array(r)
    n_vs_t[t] = np.array(n)
    F_vs_t[t] = np.array(F_tot)

# write data to disk
pickle.dump( [parameters, r_vs_t,F_vs_t,n_vs_t], open( outfile, "wb" ) )
print("--- %s seconds ---" % (time.time() - start_time))