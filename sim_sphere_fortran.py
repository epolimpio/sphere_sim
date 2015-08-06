import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri
from scipy.spatial import Delaunay
#import cPickle as pickle # For Python 3 is just pickle
import pickle
from myutils import readConfigFile, calcRotationMatrix
from fortran_lib import simforces, stripack
from datetime import datetime
import time

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

def performRotation(axis, r, n):
    """
    Perform the rotation of the vectors r and n around axis,
    which is defined, in the main program, by the angular momentum
    This warraties that the rotation is around z-axis
    """

    # unitary vector of the axis
    L = np.sqrt(np.sum(axis**2,0))
    axis = axis/L

    # z direction
    e_z = np.array([0,0,1])

    # fisrt we define the angle of rotation
    angle = np.arccos(axis.dot(e_z))

    # then we define the rotation matrix
    rot_axis = np.cross(axis,e_z)
    R = calcRotationMatrix(rot_axis, angle)
    
    # restart e_r
    N = r.shape[1]

    for i in range(0,N):
        r[:,i] = R.dot(r[:,i])
        n[:,i] = R.dot(n[:,i])
        
    e_r=get_e_r(r)

    return r, n, e_r

# ----- MAIN FUNCTION ----- #

def main(parameters):
    """
    Execute the code based on the parameters provided
    """

    # ----- PARAMETERS ----- #
    # time counter
    start_time = time.time()

    # number of bins for polar graphs
    n_bins_polar = 40

    # number of particles
    N = int(parameters['N'])

    # simulation parameters
    nu_0 = float(parameters['nu_0'])
    J = float(parameters['J'])
    eta_n = float(parameters['eta_n'])
    n_steps = int(parameters['n_steps'])
    n_save = int(parameters['n_save'])
    ftype = int(parameters['ftype'])

    # timestep
    dt = float(parameters['dt'])

    # outfile for data
    tnow = datetime.now()
    outfile_video = datetime.strftime(tnow,str(parameters['outfile_video']))
    outfile_analysis = datetime.strftime(tnow,str(parameters['outfile_analysis']))
    outfile_postrotation = datetime.strftime(tnow,str(parameters['outfile_postrotation']))

    # calculate radius of sphere from the packing parameter
    phi_pack = float(parameters['phi_pack'])
    rho=np.sqrt(N/4/phi_pack) # sphere radius

    # ----- PARAMETERS TO CONTROL ROTATION ----- #

    # Internal parameters that define when the system is in the stationary state
    rotated = False
    # number of times below the value to trigger rotation
    max_cnt_thr = int(0.2 // dt)
    # number of frames without being below the threshold to reset counter
    max_interval = int(0.1 // dt)
    # threshold of the norm of the change in angular momentum between frames
    # counting starts for variation below this value
    thr_ps_change = 5e-3*dt
    # start counters
    count_threshold = 0
    count_interval = 0

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

    # parameters data
    p_angular=[np.zeros(3),]*n_steps
    av_pairs_dist = [0.0,]*n_steps

    # relax for 2*tau
    relax_time = int(1 // dt) + 1

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

        # caculate total forces (done in FORTRAN)
        if (t<0) or (ftype == 0):
            F_tot, stress=simforces.calc_force_elastic(r,n,nu_0)
        else:
            F_tot, stress=simforces.calc_force_hooke(r,n,nu_0, list_+1)
        
        # save data before integration
        if rotated:
            # get unit vectors e_phi and e_th
            e_phi = get_e_phi(r)
            e_th = get_e_th(r)

            # save the data before integration
            stress_vs_t[t_after_rotation],bins = makePolarData(r, stress, n_bins_polar)
            
            # get force in spherical coordinates
            force = np.zeros((3,N))
            for i in range(0,N):
                force[0,i] = F_tot[:,i].dot(e_r[:,i]) # F_r
                force[1,i] = F_tot[:,i].dot(e_th[:,i]) # F_th
                force[2,i] = F_tot[:,i].dot(e_phi[:,i]) # F_phi
            force_vs_t[t_after_rotation],bins = makePolarData(r, force, n_bins_polar)

            # get directions in spherical coordinates
            direction = np.zeros((3,N))
            for i in range(0,N):
                direction[0,i] = n[:,i].dot(e_r[:,i]) # n_r
                direction[1,i] = n[:,i].dot(e_th[:,i]) # n_th
                direction[2,i] = n[:,i].dot(e_phi[:,i]) # n_phi
            direction_vs_t[t_after_rotation],bins = makePolarData(r, direction, n_bins_polar)            

            # save the number of particles per bin
            density_vs_t[t_after_rotation], bins = makePolarData(r, np.ones((1,N)), n_bins_polar, 'sum')

            t_after_rotation += 1
        
        # integrate equations of motion for each particle
        for i in range(0,N):
            # Calculate force in the local plane at position <r_i>
            F_tot_plane[:,i]=F_tot[:,i]-np.inner(F_tot[:,i],e_r[:,i])*e_r[:,i]
            
            # Calculate unit vector of total force
            f=F_tot_plane[:,i]/np.sqrt(np.sum(F_tot_plane[:,i]**2))
            
            # Calculate rotation interaction
            e_rot=np.cross(n[:,i],f)
            
            # Calculate new direction n
            dn=dt*(J*np.inner(e_r[:,i], e_rot) + eta_n*np.random.randn())*np.cross(e_r[:,i],n[:,i])
            n[:,i] = n[:,i] + dn
            n[:,i]=n[:,i]/np.sqrt(np.sum(n[:,i]**2))

            # Add to angular momentum
            ps += np.cross(r[:,i],F_tot_plane[:,i])

            # Calculate new r
            dr = dt*F_tot_plane[:,i]
            r[:,i] = r[:,i] + dr        

        # calculate unit vector <n_i> in the plane at the new position <r_i(t+dt)>
        r = constrain_to_sphere(r, rho)
        e_r=get_e_r(r)
        for i in range(0,N):
            n[:,i] = n[:,i] - np.inner(n[:,i],e_r[:,i])*e_r[:,i]
            n[:,i]=n[:,i]/np.sqrt(np.sum(n[:,i]**2))

        # Calculate neighbors
        list_ = getDelaunayTrianglesOnSphere(r)
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


        # Check if is rotating to perform global rotation
        if t>relax_time and not rotated:
            # rolling average of the angular momentum vector
            ps_av = np.average(np.array(p_angular[(t-relax_time):t]),0)
            # difference with the previous average (increment)
            dif_ps = ps_av - ps_comp

            # Check if the normalized size of the increment is below threshold
            if np.sqrt(np.sum(dif_ps**2))/N/rho/nu_0 < thr_ps_change:

                # reset interval counting 
                # (to check how long it will be needed for the next)
                count_interval = 0
                
                # If it achieved the threshold counting, rotate, else keep counting
                if count_threshold > max_cnt_thr:
                    print("Rotating at time %d"%(t+1))
                    # rotate the relevant vectors
                    r, n, e_r = performRotation(ps_av, r, n)
                    
                    # start variable after rotation
                    rotation_time = t
                    rot_axis = ps_av
                    t_after_rotation = 0
                    
                    # start variables to be saved
                    t_until_end = n_steps + relax_time - t
                    stress_vs_t = [np.zeros((9, n_bins_polar)),]*t_until_end
                    force_vs_t = [np.zeros((3, n_bins_polar)),]*t_until_end
                    direction_vs_t = [np.zeros((3, n_bins_polar)),]*t_until_end
                    density_vs_t = [np.zeros(n_bins_polar),]*t_until_end

                    rotated = True
                else:
                    count_threshold += 1
            else:
                # increase the interval to wait for the next crossing
                count_interval += 1
                if count_interval > max_interval:
                    # too many steps with no crossing, reset counting
                    count_threshold = 0

            ps_comp = ps_av
        elif t == relax_time:
            # For the first comparison
            ps_comp = np.average(np.array(p_angular),0)

    # write data to disk
    pickle.dump( [parameters, r_vs_t,F_vs_t,n_vs_t,rotated], open( outfile_video, "wb" ) )
    pickle.dump( [parameters, p_angular, av_pairs_dist,rotated], open( outfile_analysis, "wb" ) )
    if rotated:
        pickle.dump( [parameters, rotation_time, rot_axis, stress_vs_t, force_vs_t, direction_vs_t, density_vs_t, bins], open( outfile_postrotation, "wb" ) )
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == "__main__":
    """
    Run the main when externally called
    """
    parameters = readConfigFile('parameters.ini')
    main(parameters)