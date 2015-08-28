
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
#import cPickle as pickle # For Python 3 is just pickle
import pickle
from datetime import datetime
from myutils import readConfigFile, sigmoidalFunction
from metadata_lib import findFilesWithParameters, transformParameters
import sys
from os.path import join
import itertools
from scipy.optimize import curve_fit
from scipy.fftpack import fft
from scipy.signal import blackman
from scipy import stats

def openPickleFile(outfile):

    try:
        data=pickle.load(open(outfile, "rb" ) )
    except:
        print("ERROR!")
        sys.exit("Could not read file. Make sure that the file %s exists." % outfile)

    return data

# metadata path
metadata = './data/metadata.dat'
path = './data'

# start parameters with parameters.ini
# there we should set the run parameters (dt, n_save, n_steps) 
parameters = readConfigFile('parameters.ini')
parameters = transformParameters(parameters)

# simulation parameters
N = parameters['N']
n_steps = parameters['n_steps']
phi_pack = float(parameters['phi_pack'])
rho=np.sqrt(N/4/phi_pack) # sphere radius

# timestep
dx = parameters['dx']

# Now we choose the parameters we want to compare in the graph
var1_str = 'nu_0' 
val1_label = r'$\nu_0$'
var1_val = 0.1*np.arange(0,11)+1
var2_str = 'J'
val2_label = r'$J$'
var2_val = 0.05*np.arange(1,11)

var1, var2 = np.meshgrid(var1_val, var2_val)

n1, n2 = var1.shape

# start data
c_time = np.zeros(var1.shape)
max_val = np.zeros(var1.shape)

# Prepare figure
fig1 = plt.figure(1)
ax1 = plt.subplot(1,1,1)
fig2 = plt.figure(2)
ax2 = plt.subplot(1,1,1)

for i in range(0,n1):
    for j in range(0,n2):
        parameters[var1_str] = var1[i,j]
        parameters[var2_str] = var2[i,j]
        files,dates = findFilesWithParameters(metadata, parameters)
        n_files = len(files)
        # print(var1[i,j], var2[i,j], n_files)
        if n_files == 0:
            print('WARNING! No files found for var1 = {0}, var2 = {1}'.format(var1[i,j], var2[i,j]))
        # Calculate the average over the files
        nu_0 = parameters['nu_0']
        dt = dx/nu_0
        t = np.arange(n_steps)*dt
        p_mean = [0.0,]*n_steps
        p_std = [0.0,]*n_steps
        for file_ in files:
            filename = join(path, file_)
            data = openPickleFile(filename)
            p_angular=data[1]
            for step in range(0,n_steps):
                ps = np.sqrt(p_angular[step].dot(p_angular[step]))/N/rho/nu_0
                p_mean[step] += ps/n_files
                p_std[step] += ps**2/n_files 

        # Calculate mean and standard deviation
        p_mean = np.array(p_mean)
        p_std = np.array(p_std)
        p_std = p_std - p_mean**2

        popt, pcov = curve_fit(sigmoidalFunction, t, p_mean)
        a, b, c = popt

        c_time[i,j] = 1/c
        max_val[i,j] = a/b

contour1 = ax1.contourf(var1, var2, c_time)
cb1 = fig1.colorbar(contour1)
fig1.suptitle('Characteristic time', size=18, weight='bold')
ax1.set_xlabel(val1_label, size=18)
ax1.set_ylabel(val2_label, size=18)
ax1.tick_params(labelsize=16)
cb1.ax.tick_params(labelsize=16)

contour2 = ax2.contourf(var1, var2, max_val)
cb2 = fig2.colorbar(contour2)
fig2.suptitle('Maximum value order parameter', size=18, weight='bold')
ax2.set_xlabel(val1_label, size=18)
ax2.set_ylabel(val2_label, size=18)
ax2.tick_params(labelsize=16)
cb2.ax.tick_params(labelsize=16)

plt.show()