#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri
from scipy.spatial import Delaunay
import pickle
from myutils import readConfigFile, calcRotationMatrix
from fortran_lib import simforces, stripack
from datetime import datetime
import time
import metadata_lib
import random
import sys

# ----- FUNCTIONS ----- #

def get_e_r(r):
    """
    Get the cartesian coordinates of the spherical coordinate
    unit vector e_r
    """    
    l=np.sqrt(np.sum(r**2,0))
    return r/l

def get_e_th(r):
    """
    Get the cartesian coordinates of the spherical coordinate
    unit vector e_theta
    Theta here is the polar angle measured from the z-axis
    """ 
    l=np.sqrt(np.sum(r**2,0))
    th=np.arccos(r[2,:]/l)
    phi=np.arctan2(r[1,:],r[0,:])
    return np.array([np.cos(th)*np.cos(phi),np.cos(th)*np.sin(phi),-np.sin(th)])
    
def get_e_phi(r):
    """
    Get the cartesian coordinates of the spherical coordinate
    unit vector e_phi
    Phi here is the azimuthal angles measured from the x-axis
    """ 

    phi=np.arctan2(r[1,:],r[0,:])
    return np.array([-np.sin(phi),np.cos(phi),np.zeros(phi.size)])

def constrain_to_sphere(r, rho):
    """
    Contrain the vector r to the sphere of radius rho
    """ 
    l=np.sqrt(np.sum(r**2,0))
    return r/l*rho

def makePolarData(r, data, n_bins, calc_type='avg'):
    """
    INPUT: r(3,N) -> array with 3 coordinates for N particles
    INPUT: data(m, N) -> value of m parameters for each particle
    INPUT: dtheta -> size of the bin
    INPUT: calc_type -> string that defines how to calculate the data
           Possible values:
           'avg' (default) -> calculate the average of the data
           'sum' -> calculate the sum of the data
           'std' -> calculate the standard deviation of the data 
    OUTPUT: polar_data(m,n_bins) -> data distributed according to the polar angle
    If no data is found for a bin the default value is zero
    """ 
    # get the number of particles
    N = r.shape[1]
    m = data.shape[0]

    # calculate bins
    bins = np.linspace(0, np.pi, n_bins+1)
    
    # calculate polar angles
    mod_r = np.sqrt(np.sum(r**2, 0))
    theta = np.arccos(r[2,:]/mod_r)

    # place the data in the bins
    polar_data = np.zeros((m,n_bins))
    for i in range(0,n_bins):
        # get the boolean array associated with the bin
        edge1 = bins[i]
        edge2 = bins[i+1]
        if (i+1) == n_bins:
            aux = np.logical_and(
                    np.greater_equal(theta,edge1),
                    np.less_equal(theta,edge2))
        else:
            aux = np.logical_and(
                    np.greater_equal(theta,edge1),
                    np.less(theta,edge2))
        if np.sum(aux.astype(np.int)) > 0:
            if calc_type.lower() == 'sum':
                polar_data[:,i] = np.sum(data[:,aux],axis=1)
            elif calc_type.lower() == 'std':
                polar_data[:,i] = np.std(data[:,aux],axis=1)
            else: # 'avg'
                polar_data[:,i] = np.average(data[:,aux],axis=1)
        else:
            polar_data[:,i] = np.nan

    return polar_data, bins

def getVoronoiOnSphere(r):
    """
    Use STRIPACK to get the Delaunay Triangles List,
    assuming that the points are on a sphere.

    INPUT: r(3,N) -> array with 3 coordinates for N particles
    OUTPUT: tri_list(3,2*N-4) 
    """
    # Get the parameters (we do not assume unit sphere)
    N = r.shape[1]
    mod_r = np.sqrt(np.sum(r**2,0))
    r = r/mod_r
    x = r[0,:]
    y = r[1,:]
    z = r[2,:]

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
            out_ltri = ltri.T - 1

            # Get Voronoi baricenters
            lptr,lnew,ltri,listc,nb,xc,yc,zc,rc,ier = stripack.crlist(N, x, y, z, list_, lend, lptr, lnew)
            if ier == 0:
                baricenters = np.array([xc, yc, zc])
                # Construct polygons dictionary with indices of baricenter
                # also outputs all the edges pairs (for connection baricenters)
                # and the areas associated to each cell
                out_polygon_dict = {}
                pairs = []
                all_areas = np.zeros(N)
                sum_=0
                for i in range(N):
                    out_polygon_dict[i] = []
                    
                    lpl = lend[i]-1
                    index = listc[lpl]-1
                    lp = lpl

                    exit_condition = True
                    while exit_condition:
                        lp = lptr[lp]-1
                        index_prev = index
                        index = listc[lp]-1
                        if index_prev < index:
                            pairs.append([index_prev, index])
                        out_polygon_dict[i].append(index_prev)
                        all_areas[i] += stripack.areas(r[:,i], baricenters[:,index_prev], baricenters[:,index])
                        if lp == lpl:
                            exit_condition = False
                    sum_ += len(out_polygon_dict[i])

                return out_ltri, baricenters, out_polygon_dict, np.array(pairs), all_areas
            else:
                return out_ltri, None, None, None, None
 
    # Case it fails
    return None, None, None, None, None

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

def calcGlobalRotationMatrix(axis, reference = np.array([0,0,1])):

    # unitary vector of the axis
    L = np.sqrt(np.sum(axis**2,0))
    axis = axis/L

    # unitary vector of reference
    L = np.sqrt(np.sum(reference**2,0))
    reference = reference/L    

    # fisrt we define the angle of rotation
    angle = np.arccos(axis.dot(reference))

    # then we define the rotation matrix
    rot_axis = np.cross(axis, reference)
    return calcRotationMatrix(rot_axis, angle)

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
    rotate_coord = parameters['rotate_coord'] == 1
    N_fix = parameters['N_fix']
    J_chemo = parameters['J_chemo']

    # timestep
    dx = parameters['dx']
    dt = dx/nu_0

    # outfile for data
    tnow = datetime.now()
    outfile_video = datetime.strftime(tnow,str(parameters['outfile_video']))
    outfile_analysis = datetime.strftime(tnow,str(parameters['outfile_analysis']))
    outfile_postrotation = datetime.strftime(tnow,str(parameters['outfile_postrotation']))

    # calculate radius of sphere from the packing parameter
    phi_pack = float(parameters['phi_pack'])
    rho=np.sqrt(N/4/phi_pack) # sphere radius

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

    # number of bins for polar graphs (we want it even)
    n_bins_polar = int(np.floor(np.pi*np.sqrt(N)/4))
    if not n_bins_polar % 2 == 0:
        n_bins_polar -= 1

    # ----- PARAMETERS TO CONTROL ROTATION ----- #

    # Internal parameters that define when the system is in the stationary state
    rotated = False
    n_before = parameters['n_before']

    # if it will not rotate then initiate variables
    if not rotate_coord:
        rotation_time = None
        rot_axis = None

    t_after_rotation = 0

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

    # calculate unit vectors <e_r,i> and <n_i>
    e_r = get_e_r(r)
    n = np.cos(nu)*get_e_th(r) + np.sin(nu)*get_e_phi(r)

    # clear variables for storing data  
    total_saved = n_steps // n_save

    # video data
    r_vs_t=[np.zeros((3,N)),]*total_saved
    F_vs_t=[np.zeros((3,N)),]*total_saved
    n_vs_t=[np.zeros((3,N)),]*total_saved

    # Initialize distributions
    stress_polar = np.zeros((9, n_bins_polar))
    force_polar = np.zeros((3, n_bins_polar))
    direction_polar = np.zeros((3, n_bins_polar))
    density_polar = np.zeros((1,n_bins_polar))

    # parameters data
    p_angular=[np.zeros(3),]*n_steps
    av_pairs_dist = [0.0,]*n_steps
    res_force = [np.zeros(3),]*n_steps

    # relax for tau
    relax_time = int(2 // dx) + 1

    # initial neighbors list
    list_ = getDelaunayTrianglesOnSphere(r)
    # pin positions
    do_not_update_pos = []

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
            n = np.cos(nu)*get_e_th(r) + np.sin(nu)*get_e_phi(r)
            fixed_list = getDelaunayTrianglesOnSphere(r)
            if N_fix > 0:
                # fix N_pin positions
                do_not_update_pos = random.sample(range(N),N_fix)                

        # caculate total forces (done in FORTRAN)
        if (ftype == 1):
            F_tot, stress=simforces.calc_force_hooke(r,n,nu_0, list_+1)
        elif (ftype == 2):
            F_tot, stress=simforces.calc_force_hooke_break(r,n,nu_0, fanisotropy, max_dist, list_+1)
        else:
            F_tot, stress=simforces.calc_force_elastic(r,n,nu_0)
        
        # save data before integration
        if t > n_before:
            # get unit vectors e_phi and e_th
            e_phi = get_e_phi(r)
            e_th = get_e_th(r)

            # save the data before integration
            stress_dist,bins = makePolarData(r, stress, n_bins_polar)
            stress_polar += stress_dist
            
            # get force in spherical coordinates
            force = np.zeros((3,N))
            for i in range(0,N):
                force[0,i] = F_tot[:,i].dot(e_r[:,i]) # F_r
                force[1,i] = F_tot[:,i].dot(e_th[:,i]) # F_th
                force[2,i] = F_tot[:,i].dot(e_phi[:,i]) # F_phi
            force_dist,bins = makePolarData(r, force, n_bins_polar)
            force_polar += force_dist

            # get directions in spherical coordinates
            direction = np.zeros((3,N))
            for i in range(0,N):
                direction[0,i] = n[:,i].dot(e_r[:,i]) # n_r
                direction[1,i] = n[:,i].dot(e_th[:,i]) # n_th
                direction[2,i] = n[:,i].dot(e_phi[:,i]) # n_phi
            direction_dist,bins = makePolarData(r, direction, n_bins_polar)
            direction_polar += direction_dist         

            # save the number of particles per bin
            density_dist, bins = makePolarData(r, np.ones((1,N)), n_bins_polar, 'sum')
            density_polar += density_dist

            t_after_rotation += 1
        
        # integrate equations of motion for each particle
        for i in range(0,N):
            # Calculate force in the local plane at position <r_i>
            F_tot_plane[:,i]=F_tot[:,i]-np.inner(F_tot[:,i],e_r[:,i])*e_r[:,i]
            
            # Calculate unit vector of total force
            f=F_tot_plane[:,i]/np.sqrt(np.sum(F_tot_plane[:,i]**2))
            
            # Calculate rotation interaction
            e_rot=np.cross(n[:,i],f)
            force_chemo = np.zeros(3)
            if t >= 0 and chemoatt:
                for point_chemo in chemo_points:
                    dir_chemo = point_chemo - r[:,i]
                    dist_to_chemo = np.sqrt(np.sum(dir_chemo**2))
                    force_chemo += np.cross(n[:,i],dir_chemo/dist_to_chemo)*np.exp(-dist_to_chemo)

            # Calculate new direction n
            dn=dt*(np.inner(e_r[:,i], J*e_rot + J_chemo*force_chemo) +
                eta_n*np.random.randn())*np.cross(e_r[:,i],n[:,i])
            n[:,i] = n[:,i] + dn
            n[:,i]=n[:,i]/np.sqrt(np.sum(n[:,i]**2))

            # Add to angular momentum
            ps += np.cross(r[:,i],F_tot_plane[:,i])

            # Calculate new r
            if not i in do_not_update_pos:
                dr = dt*F_tot_plane[:,i]
                r[:,i] = r[:,i] + dr        

        # calculate unit vector <n_i> in the plane at the new position <r_i(t+dt)>
        r = constrain_to_sphere(r, rho)
        e_r=get_e_r(r)
        for i in range(0,N):
            n[:,i] = n[:,i] - np.inner(n[:,i],e_r[:,i])*e_r[:,i]
            n[:,i]=n[:,i]/np.sqrt(np.sum(n[:,i]**2))

        # Calculate neighbors
        if (t<0) or update_nn:
            list_ = getDelaunayTrianglesOnSphere(r)
        else:
            list_ = fixed_list
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

        # Rotate at the determined step
        if rotate_coord:
            if t==n_before and not rotated:
                # we rotate around the average angular momentum (averaged for time tau)
                ps_av = np.average(np.array(p_angular[(t-relax_time):t]),0)
                print("Rotating at step %d"%(t+1))
                # rotate the relevant vectors
                R = calcGlobalRotationMatrix(ps_av)
                for i in range(0,N):
                    r[:,i] = R.dot(r[:,i])
                    n[:,i] = R.dot(n[:,i])
                
                e_r = get_e_r(r)

                if chemoatt:
                    for point_chemo in chemo_points:
                        point_chemo = R.dot(point_chemo)

                # start variable after rotation
                rotation_time = t
                rot_axis = ps_av
                rotated = True

    # write data to disk
    pickle.dump( [parameters, r_vs_t,F_vs_t,n_vs_t,rotated,pairs,do_not_update_pos], open( outfile_video, "wb" ) )
    pickle.dump( [parameters, p_angular, av_pairs_dist, res_force, chemo_points, rotated], open( outfile_analysis, "wb" ) )

    if rotated or not rotate_coord:
        pickle.dump( [parameters, rotation_time, rot_axis,
            stress_polar/t_after_rotation, force_polar/t_after_rotation,
            direction_polar/t_after_rotation, density_polar/t_after_rotation, bins],
            open( outfile_postrotation, "wb" ) )
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == "__main__":
    """
    Run the main when externally called
    """
    parameters = readConfigFile('parameters.ini')
    main(parameters)
    print('Updating metadata...')
    metadata_lib.generateMetadata('./data')
    print('Done')