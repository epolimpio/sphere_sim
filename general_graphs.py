import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
#import cPickle as pickle # For Python 3 is just pickle
import pickle
from datetime import datetime
from myutils import readConfigFile
from metadata_lib import findFilesWithParameters, transformParameters
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
    if files:
        data = openPickleFile(join('./data', files[-1]))
        dtime = dates[-1]
    else:
        sys.exit("No file found with the parameters in parameters.ini")

print(dtime)

# Get the data
parameters=transformParameters(data[0])
p_angular=data[1]
av_dist_pairs=data[2]
res_force = np.array(data[3])
chemo_points = data[4]
print(chemo_points)
rotated = data[5]

# parameters
N = parameters['N']
n_steps = parameters['n_steps']
n_save = parameters['n_save']
phi_pack = parameters['phi_pack']
rho=np.sqrt(N/4/phi_pack) # sphere radius
nu_0 = parameters['nu_0']
dt = parameters['dx']/nu_0
rotate_coord = parameters['rotate_coord'] == 1

# If rotated, plot rotated data
if rotated or not rotate_coord:
    rotated_filename = datetime.strftime(dtime, str(parameters['outfile_postrotation']))
    data = openPickleFile(rotated_filename)
    rotation_time = data[1]
    rot_axis = data[2]
    stress_polar = data[3]
    force_polar = data[4]
    direction_polar = data[5]
    density_polar = data[6]

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

# Plot Resulting force
fig_num += 1
fig = plt.figure(fig_num)
ax = plt.subplot(111)

ax.plot(t, res_force)
ax.legend([r'$F_x$', r'$F_y$', r'$F_z$'])
plt.xlabel(r'Time ($\tau$)')
plt.ylabel(r'Resulting velocity ($\sigma \tau^{-1}$)')

if rotated or not rotate_coord:
    # Plot average stress
    fig_num += 1
    fig = plt.figure(fig_num)
    ax = plt.subplot(111)

    n_bins = stress_polar.shape[1]
    bin_edges = np.linspace(0,np.pi,n_bins+1)
    theta = 0.5*(bin_edges[:-1] + bin_edges[1:])
    ax.plot(theta, stress_polar[0,:], label=r'$\Sigma_{rr}$')
    ax.plot(theta, stress_polar[1,:], label=r'$\Sigma_{\theta r}$')
    ax.plot(theta, stress_polar[2,:], label=r'$\Sigma_{\phi r}$')
    ax.plot(theta, stress_polar[4,:], label=r'$\Sigma_{\theta \theta}$')
    ax.plot(theta, stress_polar[5,:], label=r'$\Sigma_{\phi \theta}$')
    ax.plot(theta, stress_polar[8,:], label=r'$\Sigma_{\phi \phi}$')
    ax.legend()
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'Average stress')

    # Plot directions
    fig_num += 1
    fig = plt.figure(fig_num)
    ax = plt.subplot(111)

    ax.plot(theta, direction_polar[0,:], label=r'$\hat{n}_r$')
    ax.plot(theta, direction_polar[1,:], label=r'$\hat{n}_{\theta}$')
    ax.plot(theta, direction_polar[2,:], label=r'$\hat{n}_{\phi}$')

    ax.legend()
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'Direction')

    # Plot forces
    fig_num += 1
    fig = plt.figure(fig_num)
    ax = plt.subplot(111)

    ax.plot(theta, force_polar[0,:], label=r'$\vec{F}_r$')
    ax.plot(theta, force_polar[1,:], label=r'$\vec{F}_{\theta}$')
    ax.plot(theta, force_polar[2,:]/rho/np.sin(theta), label=r'$\vec{F}_{\phi}$')

    ax.legend()
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'Force')

    # Plot Density
    fig_num += 1
    fig = plt.figure(fig_num)
    ax = plt.subplot(111)

    density_polar = np.reshape(density_polar, n_bins)
    density_polar[np.isnan(density_polar)] = 0
    ax.plot(theta, density_polar/(2*np.pi*rho**2*np.sin(theta)*abs(bin_edges[0]-bin_edges[1])))

    plt.xlabel(r'$\theta$')
    plt.ylabel(r'Density')

plt.show()