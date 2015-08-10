import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import cPickle as pickle

# this is a test, to see whether I can get the same results using a unit vector
# as measure of direction rather than an angle

R=0.4
F_th_0=0.2
K=2
L=4
dt=0.1
mu_th=0.5
mu_F=1.0
eta_n=0.5


#R=0.6
#N=16

R=0.4
N=36

outfile='test_vector1.p'

# assign random x,y positions, while z=0
r=L*np.random.rand(N,3)
r[:,2]=0
# assign random direction unit vectors n
th=2*np.pi*np.random.rand(N)

#r=np.array([[1.5,0.5,0],[0.5,1.5,0]])
#th=np.array([np.pi/2,0])

n=np.zeros((N,3))
for i in range(0,N):
    n[i,:]=np.array([np.cos(th[i]),np.sin(th[i]),0.])

r_vs_t=[]
F_vs_t=[]
n_vs_t=[]

e_r=np.array([0,0,1])

for t in range(0,300):
    print t
    
    F_r=np.zeros((N,3))
    F_th=np.zeros((N,3))
    F_tot=np.zeros((N,3))
    F_tot_plane=np.zeros((N,3))
    
    for i in range(0,N):
        F_th[i,:]=F_th_0*n[i,:]
        for j in range(i+1,N):
            # now calculate smallest distance between particle j and all (mirror) particles i
            R_ij_min=L
            for l in [-1,0,1]:
                for k in [-1,0,1]:
                    X0=l*L
                    Y0=k*L
                    rr=np.array([l*L,k*L,0.])+r[i,:]
                    r_ij=rr-r[j,:]
                    R_ij=np.sqrt(np.sum(r_ij**2))
                    if R_ij<R_ij_min:
                        R_ij_min=R_ij
                        r_ij_min=r_ij
            r_ij=r_ij_min
            R_ij=R_ij_min
            # check if particles i and j overlap
            f=2*R-R_ij_min
            if f>0:
                # if so, calculate repulsive force on each particle
                F_ij=K*r_ij_min/R_ij_min*f
                F_ji=-F_ij
            else:
                F_ij=np.array([0,0,0])
                F_ji=np.array([0,0,0])
            # add to total repulsive force    
            F_r[i,:]=F_r[i,:]+F_ij
            F_r[j,:]=F_r[j,:]+F_ji

    # now calculate the total force on each particle and total force in plane  
    F_tot=F_r+F_th
    for i in range(0,N):
        F_tot_plane[i,:]=F_tot[i,:]-np.inner(F_tot[i,:],e_r)*r[i,:]

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

    r=r+dt*mu_F*F_tot_plane

    r_vs_t.append(np.array(r))
    n_vs_t.append(np.array(n))
    F_vs_t.append(np.array(F_tot))
    
    for i in range(0,N):
        for j in range(0,2):
            if r[i,j]<0:
                r[i,j]=r[i,j]+L
            elif r[i,j]>L:
                r[i,j]=r[i,j]-L

pickle.dump( [r_vs_t,F_vs_t,n_vs_t], open( outfile, "wb" ) )
