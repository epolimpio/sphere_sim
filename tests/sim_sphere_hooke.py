import numpy as np
import matplotlib.pyplot as plt
#import cPickle as pickle # For Python 3 is just pickle
import pickle

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

def start_histogram(R, dr):

    num_bins = np.floor(R/dr)+1
    return np.zeros(num_bins)

def place_on_histogram(r, dr, histogram):

    bin = np.floor(r/dr)
    if bin < histogram.size:
        histogram[bin] += 1

    return histogram

# ----- PARAMETERS ----- #
# number of particles
N=100

# simulation parameters
F_th_0=0.2 #
mu_th=1.0 #
mu_F=1.0 #
eta_n=0.5 #
K=2 # elastic constant of sphere repulsion
R=0.4 # radius of spheres

# timestep
dt=0.1

# outfile for data
outfile='sphere_data_hooke.p'

# calculate radius of sphere... kind of ugly due to random parameter <L>
L=0.4 # "typical length" of the particle 
rho=np.sqrt(N*L**2/(4*np.pi)) # sphere radius
dr_hist = 0.04
histogram = start_histogram(rho, dr_hist)

# ----- INITIALIZE RANDOM POSITION AND ORIENTATION ----- #
# get random positions (i.e. <th>,<phi> in spherical coordinates) for all particles
th=np.pi*np.random.rand(N)
phi=2*np.pi*np.random.rand(N)
# get random angle for unit vector <n> in the local plane
nu=2*np.pi*np.random.rand(N)

r=np.zeros((N,3))
n=np.zeros((N,3))
e_r=np.zeros((N,3))

# calculate positions <r_i> and unit vectors <e_r,i> and <n_i>
for i in range(0,N):
    r[i,:]=rho*np.array([np.sin(th[i])*np.cos(phi[i]),np.sin(th[i])*np.sin(phi[i]),np.cos(th[i])])
    e_r[i,:]=get_e_r(r[i,:])
    n[i,:]=np.cos(nu[i])*get_e_th(r[i,:]) + np.sin(nu[i])*get_e_phi(r[i,:])

n_steps = 200
# clear variables for storing data    
r_vs_t=[np.array([]),]*n_steps
F_vs_t=[np.array([]),]*n_steps
n_vs_t=[np.array([]),]*n_steps
number_of_neighbors = [0.0,]*n_steps
number_of_new_neighbors = [0.0,]*n_steps

old_list_of_neighbors = {}
for i in range(0,N):
    old_list_of_neighbors[i] = []

# ----- RUN FOR SEVERAL TIMES ----- #
t_stab = 20
for t in range(0,n_steps):
    print(t)

    # start with zeros
    F_th=np.zeros((N,3))
    F_r=np.zeros((N,3))
    F_tot=np.zeros((N,3))
    F_tot_plane=np.zeros((N,3))
    NN=np.zeros(N)
    NN_new=np.zeros(N)
    list_of_neighbors = {}
    for i in range(0,N):
        list_of_neighbors[i] = []

    # calculate forces
    for i in range(0,N):
        # calculate active force in direction <n_i>
        F_th[i,:]=F_th_0*n[i,:]
        for j in range(i+1,N):
            # now calculate distance <R> between particles i and j
            r_ij=r[i,:]-r[j,:]
            R_ij=np.sqrt(np.sum(r_ij**2))
            # check if particles i and j overlap
            f=2*R-R_ij
            if f>=-R/2:
                # if so, calculate repulsive force on each particle
                F_ij=K*r_ij/R_ij*f
                F_ji=-F_ij
                
                # neighbors
                NN[i]+=1
                NN[j]+=1
                list_of_neighbors[i].append(j)
                list_of_neighbors[j].append(i)
                if not j in old_list_of_neighbors[i]:
                    NN_new[i] += 1
                    NN_new[j] += 1
            else:
                F_ij=np.array([0,0,0])
                F_ji=np.array([0,0,0])
            # add to total repulsive force    
            F_r[i,:]=F_r[i,:]+F_ij
            F_r[j,:]=F_r[j,:]+F_ji
            if t >= t_stab:
                histogram = place_on_histogram(R_ij, dr_hist, histogram)

    # caculate total of active and repulsive forces
    F_tot=F_r+F_th
    # calculate force in the local plane at position <r_i>
    for i in range(0,N):    
        F_tot_plane[i,:]=F_tot[i,:]-np.inner(F_tot[i,:],e_r[i,:])*e_r[i,:]
    
    # integrate equation of motion for direction n_i
    for i in range(0,N):
        # first calculate unit vector of total force
        f=F_tot_plane[i,:]/np.sqrt(np.sum(F_tot_plane[i,:]**2))
        # next, calculate axis of rotation
        e_rot=np.cross(n[i,:],f)
        e_rot_L=np.sqrt(np.sum(e_rot**2))
        if e_rot_L>1e-4:
            e_rot=e_rot/e_rot_L
        # then, calculate absolute angle A between n and f
        B=np.inner(n[i,:],f)
        if B>1.0:
            A=0.0
        elif B<-1.0:
            A=np.pi
        else:             
            A=np.arccos(np.inner(n[i,:],f))
        # finally, calculate new direction n
        dn=dt*(mu_th*A + eta_n*np.random.randn())*np.cross(e_rot,n[i,:])
        n[i,:] = n[i,:] + dn
        n[i,:]=n[i,:]/np.sqrt(np.sum(n[i,:]**2))

    # integrate equation of motion for position <r_i>
    for i in range(0,N):  
        dr = dt*mu_F*F_tot_plane[i,:]
        r[i,:] = r[i,:] + dr        
        r[i,:] = constrain_to_sphere(r[i,:])

    # calculate unit vector <n_i> in the plane at the new position <r_i(t+dt)>
    for i in range(0,N):
        e_r[i,:]=get_e_r(r[i,:])
        n[i,:] = n[i,:] - np.inner(n[i,:],e_r[i,:])*e_r[i,:]
        n[i,:]=n[i,:]/np.sqrt(np.sum(n[i,:]**2))

    # store data for <r_i>, <N-i> and <F_tot,i>
    r_vs_t[t] = np.array(r)
    n_vs_t[t] = np.array(n)
    F_vs_t[t] = np.array(F_tot)
    number_of_neighbors[t] = np.average(NN)
    number_of_new_neighbors[t] = np.average(NN_new)

    old_list_of_neighbors = list_of_neighbors

print(histogram/(n_steps-t_stab))

# write data to disk
pickle.dump( [r_vs_t,F_vs_t,n_vs_t], open( outfile, "wb" ) )
pickle.dump( [number_of_neighbors, number_of_new_neighbors, histogram/(n_steps-t_stab)], open( 'neighbors_hooke.p', "wb" ) )
