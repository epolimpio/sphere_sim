
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

# Now we choose the parameters we want to compare in the graph
# conditions: N, phi, eta, nu, J
all_N = [100]
all_phi = [1]
all_nu = [1]
all_J = [0.15, 0.6, 1.05, 1.5]
all_eta = [0.6]
all_anisotropy = [1]
all_max_dist = [0]
all_update_nn = [1]

all_comb = [all_N, all_phi, all_nu, all_J, all_eta, all_anisotropy, all_max_dist, all_update_nn]
conditions = list(itertools.product(*all_comb))

# Here we define the legend style
variables_legend = ['J']
legend_format = r'$J = {0}$'
# define the parameters of the graphs
combinations = {}
cnt = 0
for condition in conditions:
    N, phi_pack, nu_0, J, eta_n, anis, dist, update_nn = condition
    # change the parameters
    parameters['N'] = N
    parameters['nu_0'] = nu_0
    parameters['J'] = J
    parameters['eta_n'] = eta_n
    parameters['phi_pack'] = phi_pack
    parameters['fanisotropy'] = anis
    parameters['max_dist'] = dist
    parameters['update_nn'] = update_nn
    combinations[cnt] = parameters.copy()
    cnt += 1

# Prepare figure
fig = plt.figure(1)
ax = plt.subplot(2,2,1)
ax2 = plt.subplot(2,2,2)
ax3 = plt.subplot(2,2,3)
ax4 = plt.subplot(2,2,4)
fig2 = plt.figure(2)
ax_sub = plt.subplot(1,1,1)
variable = []
chract_time = []
max_val = []
for cnt in combinations:
    combination = combinations[cnt]
    files,dates = findFilesWithParameters(metadata, combination)
    n_files = len(files)

    combination = transformParameters(combination)
    # parameters
    N = combination['N']
    n_steps = combination['n_steps']
    n_save = combination['n_save']
    phi_pack = combination['phi_pack']
    rho=np.sqrt(N/4/phi_pack) # sphere radius
    nu_0 = combination['nu_0']
    eta_n = combination['eta_n']
    J = combination['J']
    dt = combination['dx']/nu_0
    t = np.arange(n_steps)*dt

    print('\nParameters:')
    print('N: {0}, phi: {1}, J: {2}, nu: {3}, eta: {4}'.format(N, phi_pack, J, nu_0, eta_n))
    print('# files: {0}'.format(n_files))

    # Calculate the average over the files
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

    # Calculate standard deviation
    p_mean = np.array(p_mean)
    p_std = np.array(p_std)
    p_std = p_std - p_mean**2

    popt, pcov = curve_fit(sigmoidalFunction, t, p_mean)
    a, b, c = popt

    variable.append(combination[variables_legend[0]])
    chract_time.append(c)
    max_val.append(a/b)

    p_mean_sub = p_mean - sigmoidalFunction(t, *popt)
    # Plot the result
    legend_data = []
    for var in variables_legend:
        legend_data.append(combination[var])
    legend_str = legend_format.format(*legend_data)

    # Plot absolute values
    line, = ax.plot(t, p_mean, label=legend_str)
    ax.plot(t, sigmoidalFunction(t, *popt), '--k')
    color = line.get_color()
    ax.fill_between(t, p_mean-p_std, p_mean+p_std, facecolor=color, alpha=0.3)

    # Plot FFT
    freq = np.fft.fftfreq(t.shape[-1])
    w = blackman(n_steps)
    sp = np.fft.fft(p_mean_sub) 
    ax3.plot(freq[0:n_steps*dt]/dt, np.abs(sp)[0:n_steps*dt])

    ax_sub.plot(t, p_mean_sub)

ax2.plot(variable, max_val)
ax4.plot(variable, chract_time)
slope, intercept, r_value, p_value, std_err = stats.linregress(np.array(variable),np.array(chract_time))
print(slope, r_value**2)


#ax.set_yscale('log')
ax.legend(loc=4, ncol=3)
ax.set_xlabel(r'Time ($\tau$)')
ax.set_ylabel(r'Order parameter')
ax2.set_xlabel(r'J ($\tau^{-1}$)')
ax2.set_ylabel(r'Order Parameter ($t \to \infty$)')
ax3.set_xlabel(r'Frequency')
ax3.set_ylabel(r'Contribution')
ax4.set_xlabel(r'J ($\tau^{-1}$)')
ax4.set_ylabel(r'Characteristic Time ($\tau$)')



plt.show()