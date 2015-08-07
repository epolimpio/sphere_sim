import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
#import cPickle as pickle # For Python 3 is just pickle
import pickle
from datetime import datetime
from myutils import readConfigFile
from metadata_lib import findFilesWithParameters
import sys
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
    
    analysis_filename = './data/sphere_data_analysis_%Y-%m-%d_%H-%M-%S.p'
    dtime = datetime(yy,mm,dd,h,m,s)
    outfile=datetime.strftime(dtime,analysis_filename)
    data = openPickleFile(outfile)
else:
    # get the last file with the parameters in parameters.ini
    parameters = readConfigFile('parameters.ini')
    metadata = './data/metadata.dat'
    files, dates = findFilesWithParameters(metadata, parameters)
    data = openPickleFile(join('./data', files[-1]))
    dtime = dates[-1]    

# Get the data
parameters=data[0]
p_angular=data[1]
av_dist_pairs=data[2]
rotated = data[3]

# parameters
N = int(parameters['N'])
n_steps = int(parameters['n_steps'])
n_save = int(parameters['n_save'])
phi_pack = float(parameters['phi_pack'])
rho=np.sqrt(N/4/phi_pack) # sphere radius
nu_0 = float(parameters['nu_0'])
dt = float(parameters['dt'])

# If rotated, plot rotated data
if rotated:
    rotated_filename = datetime.strftime(dtime, str(parameters['outfile_postrotation']))
    data = openPickleFile(rotated_filename)
    rotation_time = data[1]
    rot_axis = data[2]
    stress_vs_t = np.array(data[3])
    force_vs_t = np.array(data[4])
    direction_vs_t = np.array(data[5])
    density_vs_t = data[6]

fig_num = 0

# Plot Order Parameter
fig_num += 1
fig = plt.figure(fig_num)
ax = plt.subplot(111)

ps = [0.0,]*n_steps
for i in range(0,n_steps):
    ps[i] = np.sqrt(p_angular[i].dot(p_angular[i]))/N/rho/nu_0

t = np.arange(n_steps)*dt
ax.plot(t, ps)
plt.xlabel(r'Time ($\tau$)')
plt.ylabel(r'Order parameter')

# Plot Average pair distance
fig_num += 1
fig = plt.figure(fig_num)
ax = plt.subplot(111)

ax.plot(t, av_dist_pairs)
plt.xlabel(r'Time ($\tau$)')
plt.ylabel(r'Average pair distance ($\sigma$)')

if rotated:
    # Plot average stress
    fig_num += 1
    fig = plt.figure(fig_num)
    ax = plt.subplot(111)

    n_bins = stress_vs_t.shape[2]
    bin_edges = np.linspace(0,np.pi,n_bins+1)
    theta = 0.5*(bin_edges[:-1] + bin_edges[1:])

    ax.plot(theta, np.nanmean(stress_vs_t[:,0,:],axis = 0), label=r'$\Sigma_{rr}$')
    ax.plot(theta, np.nanmean(stress_vs_t[:,1,:],axis = 0), label=r'$\Sigma_{\theta r}$')
    ax.plot(theta, np.nanmean(stress_vs_t[:,2,:],axis = 0), label=r'$\Sigma_{\phi r}$')
    ax.plot(theta, np.nanmean(stress_vs_t[:,4,:],axis = 0), label=r'$\Sigma_{\theta \theta}$')
    ax.plot(theta, np.nanmean(stress_vs_t[:,5,:],axis = 0), label=r'$\Sigma_{\phi \theta}$')
    ax.plot(theta, np.nanmean(stress_vs_t[:,8,:],axis = 0), label=r'$\Sigma_{\phi \phi}$')
    ax.legend()
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'Average stress')

    # Plot directions
    fig_num += 1
    fig = plt.figure(fig_num)
    ax = plt.subplot(111)

    ax.plot(theta, np.nanmean(direction_vs_t[:,0,:],axis = 0), label=r'$\hat{n}_r$')
    ax.plot(theta, np.nanmean(direction_vs_t[:,1,:],axis = 0), label=r'$\hat{n}_{\theta}$')
    ax.plot(theta, np.nanmean(direction_vs_t[:,2,:],axis = 0), label=r'$\hat{n}_{\phi}$')

    ax.legend()
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'Direction')

    # Plot forces
    fig_num += 1
    fig = plt.figure(fig_num)
    ax = plt.subplot(111)

    ax.plot(theta, np.nanmean(force_vs_t[:,0,:],axis = 0), label=r'$\vec{F}_r$')
    ax.plot(theta, np.nanmean(force_vs_t[:,1,:],axis = 0), label=r'$\vec{F}_{\theta}$')
    ax.plot(theta, np.nanmean(force_vs_t[:,2,:],axis = 0), label=r'$\vec{F}_{\phi}$')

    ax.legend()
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'Force')

    # Plot Density
    fig_num += 1
    fig = plt.figure(fig_num)
    ax = plt.subplot(111)

    dens_dist = np.zeros(n_bins)
    for dist in density_vs_t:
        dist = dist.reshape(n_bins)
        dist[np.isnan(dist)] = 0
        dens_dist += dist

    ax.plot(theta, dens_dist/len(density_vs_t)/(2*np.pi*rho**2*np.sin(theta)*abs(bin_edges[0]-bin_edges[1])))

    plt.xlabel(r'$\theta$')
    plt.ylabel(r'Density')

plt.show()