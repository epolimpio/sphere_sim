#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri
import matplotlib.patches as patches
from scipy.spatial import Delaunay
import pickle
from myutils import readConfigFile, calcRotationMatrix
from fortran_lib import simforces
from scipy.spatial import Delaunay, ConvexHull
from datetime import datetime
import time
import metadata_lib
import random
import sys

# ----- FUNCTIONS ----- #

def getDelaunayTrianglesOnPlane(r):
    """
    Get Delaunay triagulation on plane

    INPUT: r(2,N) -> array with 2 coordinates for N particles on plane
    OUTPUT: tri_list(3,2*N-4) with neighbors list index
            boundary(2*N-4) boolean array with elements in the boundary
    """
    # Get the number of particles
    N = r.shape[1]
   
    # Calculate the Delaunay triangulation
    tri = Delaunay(r.T)
    tri_list = tri.simplices

    # Get all triangle in boundary
    where_boundary = np.equal(tri.neighbors, -1)
    # make an array in which the lines are true if triangle is in boundary
    tri_in_boundary = np.tile(np.greater(np.sum(where_boundary,axis=1),0),(3,1)).T

    # boundary is true for the vertices in boundary
    boundary = np.logical_and(np.logical_not(where_boundary), tri_in_boundary)

    # Indices of the points in boundary
    points_in_boundary = np.unique(tri.simplices[boundary].flatten())

    # return boolean array that is true for boundary
    return tri_list, np.in1d(np.arange(N), points_in_boundary)

def startPositions(N):

    Nx = int(np.sqrt(N))+2
    Ny = int(np.sqrt(N))+2
    N = Nx*Ny
    r = np.zeros((3,N))
    r[0,0:Nx] = 2*np.arange(Nx)
    disloc = np.sqrt(3)
    boundary = []
    boundary.extend(list(range(0,Nx)))
    for i in range(1,Ny):
        r[0,i*Nx:(i+1)*Nx] = r[0,(i-1)*Nx:i*Nx] + (2*(i % 2) - 1)
        r[1,i*Nx:(i+1)*Nx] = r[1,(i-1)*Nx:i*Nx] + disloc
        boundary.extend([i*Nx, (i+1)*Nx-1])
        if i == Ny-1:
            boundary.extend(list(range(i*Nx,(i+1)*Nx)))

    return r, Nx, Ny, boundary

# ----- MAIN FUNCTION ----- #

def main(parameters):
    """
    Execute the code based on the parameters provided
    """

    # ----- PARAMETERS ----- #
    # time counter
    start_time = time.time()

    # transform the parameters properly
    parameters = metadata_lib.transformParameters(parameters)

    # number of particles
    N = parameters['N']

    # simulation parameters
    nu_0 = parameters['nu_0']
    J = parameters['J']
    eta_n = parameters['eta_n']
    n_steps = parameters['n_steps']
    n_save = parameters['n_save']
    ftype = parameters['ftype']
    fanisotropy = parameters['fanisotropy']
    max_dist = parameters['max_dist']
    update_nn = parameters['update_nn'] == 1
    chemoatt = parameters['chemoatt'] == 1
    N_fix = parameters['N_fix']
    J_chemo = parameters['J_chemo']

    # timestep
    dx = parameters['dx']
    dt = dx/nu_0

    # outfile for data
    tnow = datetime.now()
    outfile_video = datetime.strftime(tnow,str(parameters['outfile_video_plane']))
    outfile_analysis = datetime.strftime(tnow,str(parameters['outfile_analysis_plane']))

    # calculate radius of sphere from the packing parameter
    phi_pack = float(parameters['phi_pack'])
    rho=np.sqrt(2*np.sqrt(3)*N/phi_pack) # Size of the square plane

    # include the points of chemo attractant
    if chemoatt:
        chemo_points_temp = eval(parameters['chemo_points'])
        chemo_points = []
        for point in chemo_points_temp:
            if not len(point) == 3:
                sys.exit('Error in the chemotaxis point')
            else:
                chemo_points.append(np.array(point))
    else:
        chemo_points = []

    # ----- PARAMETERS TO CONTROL ROTATION ----- #

    # Internal parameters that define when the system is in the stationary state
    rotated = False
    n_before = parameters['n_before']

    t_after_rotation = 0

    # ----- INITIALIZE RANDOM POSITION AND ORIENTATION ----- #

    # Start position in hexagonal array
    r, Nx, Ny, b = startPositions(N)
    N = Nx*Ny
    # pin positions of boundary particles
    do_not_update_pos = b
    
    # get random angle for unit vector <n> in the local plane
    nu=2*np.pi*np.random.rand(N)

    # Unit vectors for all the particles
    e_x = np.tile(np.array([1,0,0]).reshape(3,1), N)
    e_y = np.tile(np.array([0,1,0]).reshape(3,1), N)
    e_z = np.tile(np.array([0,0,1]).reshape(3,1), N)

    # calculate unit vectors
    n = np.cos(nu)*e_x + np.sin(nu)*e_y

    # clear variables for storing data  
    total_saved = n_steps // n_save

    # video data
    r_vs_t=[np.zeros((3,N)),]*total_saved
    F_vs_t=[np.zeros((3,N)),]*total_saved
    n_vs_t=[np.zeros((3,N)),]*total_saved

    # parameters data
    p_angular=[np.zeros(3),]*n_steps
    av_pairs_dist = [0.0,]*n_steps
    res_force = [np.zeros(3),]*n_steps

    # relax for tau
    relax_time = int(2 // dx) + 1

    # initial neighbors list
    list_, boundary = getDelaunayTrianglesOnPlane(r[0:2,:])
    boundary = np.in1d(np.arange(N), b) 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.add_patch(
        patches.Rectangle((0, 0), 2*(Nx-1), np.sqrt(3)*(Ny-1), fill=False)
    )
    ax.triplot(r[0,:], r[1,:], list_)
    ax.plot(r[0,:], r[1,:], 'o')
    ax.plot(r[0,boundary], r[1,boundary], 'xk')
    plt.show()

    #sys.exit('Ok')

    # ----- SIMULATION ----- #
    for step in range(0,n_steps+relax_time):

        # run time variable
        t = step-relax_time

        # Some prints to see program running
        if (t>=0) and ((t+1) % 50 == 0):
            print('Time: {t}, Order parameter: {ps}'.format(t=int(t+1), ps = np.sqrt(ps.dot(ps))/N/nu_0/rho))
        elif (t<0) and ((step+1) % 10 == 0):
            print('Relaxation Time: %d' % (step+1))

        # Restart variables
        stress = np.zeros((9,N))
        F_tot_plane=np.zeros((3,N))
        ps = np.zeros(3)

        # Restart directions after relaxation
        if t == 0:
            print("Restart directions")
            nu=2*np.pi*np.random.rand(N)
            n = np.cos(nu)*e_x + np.sin(nu)*e_y
            fig = plt.figure()
            ax = fig.add_subplot(111, aspect='equal')
            ax.add_patch(
                patches.Rectangle((0, 0), 2*(Nx-1), np.sqrt(3)*(Ny-1), fill=False)
            )
            ax.triplot(r[0,:], r[1,:], list_)
            ax.plot(r[0,:], r[1,:], 'o')
            ax.plot(r[0,boundary], r[1,boundary], 'xk')
            plt.show()

            fixed_list = list_
            fixed_boundary = boundary

            if N_fix > 0:
                # fix N_pin positions
                do_not_update_pos = random.sample(range(N),N_fix)                

        # caculate total forces (done in FORTRAN)
        boundary_int = boundary.astype(np.int)
        if (ftype == 1):
            F_tot, stress=simforces.calc_force_hooke(r, n, nu_0, list_+1)
        elif (ftype == 2):
            F_tot, stress=simforces.calc_force_hooke_break_plane(r, n, nu_0, fanisotropy, max_dist, rho, boundary_int, list_+1)
        else:
            F_tot, stress=simforces.calc_force_elastic_plane(r, n, nu_0, rho, boundary_int)
        
        F_tot_plane = np.copy(F_tot)
        F_tot_plane[2,:] = 0     
        # integrate equations of motion for each particle
        for i in range(0,N):
            
            # Calculate unit vector of total force
            f=F_tot_plane[:,i]/np.sqrt(np.sum(F_tot_plane[:,i]**2))
            
            # Calculate rotation interaction
            e_rot=np.cross(n[:,i],f)
            force_chemo = np.zeros(3)
            if t >= 0 and chemoatt:
                for point_chemo in chemo_points:
                    dir_chemo = point_chemo - r[:,i]
                    dist_to_chemo = np.sqrt(np.sum(dir_chemo**2))
                    force_chemo += np.cross(n[:,i],dir_chemo/dist_to_chemo)*np.exp(-dist_to_chemo/2)

            # Calculate new direction n
            dn=dt*(np.inner(e_z[:,i], J*e_rot + J_chemo*force_chemo) +
                eta_n*np.random.randn())*np.cross(e_z[:,i],n[:,i])
            n[:,i] = n[:,i] + dn
            n[:,i]=n[:,i]/np.sqrt(np.sum(n[:,i]**2))

            # Add to angular momentum
            ps += np.cross(r[:,i],F_tot_plane[:,i])

            # Calculate new r
            if not i in do_not_update_pos:
                dr = dt*F_tot_plane[:,i]
                r[:,i] = r[:,i] + dr        

        # calculate unit vector <n_i> in the plane at the new position <r_i(t+dt)>
        n[2,:] = 0
        n = n/np.sqrt(np.sum(n**2,0))

        #Calculate neighbors
        # if (t<0) or update_nn:
        #     list_, boundary1 = getDelaunayTrianglesOnPlane(r[0:2,:])
        # else:
        #     list_ = fixed_list
        #     boundary = fixed_boundary
        
        pairs_dist, pairs = simforces.calc_pairs_dist(r, list_+1)

        # Store data
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
            res_force[t] = np.sum(F_tot, axis=1)

    # write data to disk
    pickle.dump( [parameters, r_vs_t,F_vs_t,n_vs_t,rotated,pairs,do_not_update_pos], open( outfile_video, "wb" ) )
    pickle.dump( [parameters, p_angular, av_pairs_dist, res_force, chemo_points, rotated], open( outfile_analysis, "wb" ) )

    print("--- %s seconds ---" % (time.time() - start_time))

    list_, boundary = getDelaunayTrianglesOnPlane(r[0:2,:])
    boundary = np.in1d(np.arange(N), b) 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.add_patch(
        patches.Rectangle((0, 0), 2*(Nx-1), np.sqrt(3)*(Ny-1), fill=False)
    )
    ax.triplot(r[0,:], r[1,:], list_)
    ax.plot(r[0,:], r[1,:], 'o')
    ax.plot(r[0,boundary], r[1,boundary], 'xk')
    plt.show()

if __name__ == "__main__":
    """
    Run the main when externally called
    """
    parameters = readConfigFile('parameters.ini')
    main(parameters)
    print('Updating metadata...')
    metadata_lib.generateMetadata('./data', plane = True)
    print('Done')