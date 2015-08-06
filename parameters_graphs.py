import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
#import cPickle as pickle # For Python 3 is just pickle
import pickle
from datetime import datetime
import sys


def openPickleFile(filename, dtime):

    outfile=datetime.strftime(dtime,filename)

    try:
        data=pickle.load(open(outfile, "rb" ) )
    except:
        print("ERROR!")
        sys.exit("Could not read file. Make sure that the file %s exists." % outfile)

    return data

# Read the date time parameters
if len(sys.argv)>=2:
    str_date = sys.argv[1]
else:
    str_date = input("Insert date, format Y,m,d,H,M,S ->")

try:
    yy,mm,dd,h,m,s = [int(x) for x in str_date.split(',')[:6]]
except:
    print("ERROR!")
    sys.exit("Invalid date. The format must be Y,m,d,H,M,S separated by comma.")

# Open analysis file
analysis_filename = './data/sphere_data_analysis_%Y-%m-%d_%H-%M-%S.p'
dtime_data = datetime(yy,mm,dd,h,m,s)
data = openPickleFile(analysis_filename, dtime_data)

# Get the data
parameters=data[0]
p_angular=data[1]
av_dist_pairs=data[2]
rotated = data[3]

# If rotated, get rotated data
if rotated:
    rotated_filename = str(parameters['outfile_postrotation'])
    data = openPickleFile(rotated_filename, dtime_data)

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
    rotated_filename = str(parameters['outfile_postrotation'])
    data = openPickleFile(rotated_filename, dtime_data)
    rotation_time = data[1]
    rot_axis = data[2]
    stress_vs_t = np.array(data[3])
    force_vs_t = data[4]
    direction_vs_t = data[5]
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

# Plot average stress
fig_num += 1
fig = plt.figure(fig_num)
ax = plt.subplot(111)

n_bins = stress_vs_t.shape[2]
bin_edges = np.linspace(0,np.pi,n_bins+1)
theta = 0.5*(bin_edges[:-1] + bin_edges[1:])

ax.plot(theta, np.nanmean(stress_vs_t[:,0,:],axis = 0))
ax.plot(theta, np.nanmean(stress_vs_t[:,1,:],axis = 0))
ax.plot(theta, np.nanmean(stress_vs_t[:,2,:],axis = 0))
ax.plot(theta, np.nanmean(stress_vs_t[:,4,:],axis = 0))
ax.plot(theta, np.nanmean(stress_vs_t[:,5,:],axis = 0))
ax.plot(theta, np.nanmean(stress_vs_t[:,8,:],axis = 0))


plt.show()