import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri
from scipy.spatial import Delaunay
#import cPickle as pickle # For Python 3 is just pickle
import pickle
from myutils import readConfigFile 
from fortran_lib import simforces, stripack
from datetime import datetime
import time
start_time = time.time()

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
eta_n = float(parameters['eta_n'])
n_steps = int(parameters['n_steps'])
n_save = int(parameters['n_save'])
ftype = int(parameters['ftype'])

# timestep
dt = float(parameters['dt'])

# outfile for data
outfile_video = datetime.strftime(datetime.now(),str(parameters['outfile_video']))
outfile_analysis = datetime.strftime(datetime.now(),str(parameters['outfile_analysis']))

# calculate radius of sphere from the packing parameter
fpacking = float(parameters['fpacking'])
rho=np.sqrt(N/4/fpacking) # sphere radius

# ----- INITIALIZE RANDOM POSITION AND ORIENTATION ----- #
# get random positions for all particles
x = 2*np.random.rand(N) - 1
y = 2*np.random.rand(N) - 1
z = 2*np.random.rand(N) - 1
# calculate module and place them on the sphere
mod_r = np.sqrt(x**2+y**2+z**2)
x = rho*x/mod_r
y = rho*y/mod_r
z = rho*z/mod_r

r = np.vstack((x,y,z))

# get random angle for unit vector <n> in the local plane
nu=2*np.pi*np.random.rand(N)
n=np.zeros((3,N))
e_r=np.zeros((3,N))

# calculate unit vectors <e_r,i> and <n_i>
for i in range(0,N):
    e_r[:,i]=get_e_r(r[:,i])
    n[:,i]=np.cos(nu[i])*get_e_th(r[:,i]) + np.sin(nu[i])*get_e_phi(r[:,i])

# clear variables for storing data  
total_saved = n_steps // n_save

# video data
r_vs_t=[np.array([]),]*total_saved
F_vs_t=[np.array([]),]*total_saved
n_vs_t=[np.array([]),]*total_saved

# parameters data
p_angular=[np.zeros(3),]*n_steps
av_pairs_dist = [0.0,]*n_steps

# relax for 2*tau
relax_time = int(1 // dt) + 1

# ----- RUN FOR SEVERAL TIMES ----- #
for k in range(0,n_steps+relax_time):

    t = k-relax_time

    if (t>=0) and ((t+1) % 10 == 0):
        print('Time: %d' % (t+1))
    elif (k+1) % 10 == 0:
        print('Relaxation Time: %d' % (k+1))

    # caculate total of active and repulsive forces
    if (t<0) or (ftype == 0):
        F_tot=simforces.calc_force_elastic(r,n,F_th_0)
    else:
        F_tot=simforces.calc_force_hooke(r,n,F_th_0, list_+1)
    
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
    ps = np.zeros(3)
    for i in range(0,N):
        ps += np.cross(r[:,i],F_tot_plane[:,i])
        dr = dt*F_tot_plane[:,i]
        r[:,i] = r[:,i] + dr        
        r[:,i] = constrain_to_sphere(r[:,i])

    # calculate unit vector <n_i> in the plane at the new position <r_i(t+dt)>
    for i in range(0,N):
        e_r[:,i]=get_e_r(r[:,i])
        n[:,i] = n[:,i] - np.inner(n[:,i],e_r[:,i])*e_r[:,i]
        n[:,i]=n[:,i]/np.sqrt(np.sum(n[:,i]**2))

    # Calculate neighbors
    list_ = getDelaunayTrianglesOnSphere(r)
    pairs_dist, pairs = simforces.calc_pairs_dist(r, list_+1)

    if t>=0:
        if (t+1) % n_save == 0:
            # store data for <r_i>, <N-i> and <F_tot,i> for image
            index = (t+1) // n_save - 1
            r_vs_t[index] = np.array(r)
            n_vs_t[index] = np.array(n)
            F_vs_t[index] = np.array(F_tot)

        # store data for analysis (do not scale with N, so save for all t)
        p_angular[t] = ps
        av_pairs_dist[t] = np.mean(pairs_dist)

# write data to disk
pickle.dump( [parameters, r_vs_t,F_vs_t,n_vs_t], open( outfile_video, "wb" ) )
pickle.dump( [parameters, p_angular, av_pairs_dist], open( outfile_analysis, "wb" ) )
print("--- %s seconds ---" % (time.time() - start_time))