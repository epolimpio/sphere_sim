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
from os.path import join

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
# conditions: ftype, N, phi, eta, nu, J
conditions = [
[1, 200, 1, 1, 1, 0.5],
[1, 200, 1, 1, 1, 1],
[1, 200, 1, 1, 1, 2],
]

variable_compared = 'J'
legend_name = r'$J$'

# define the parameters of the graphs
combinations = {}
cnt = 0
for condition in conditions:
    [ftype, N, phi, eta, nu, J] = condition
    parameters['N'] = N
    parameters['phi_pack'] = phi
    parameters['nu_0'] = nu
    parameters['J'] = J
    parameters['eta_n'] = eta
    parameters['ftype'] = ftype
    combinations[cnt] = parameters.copy()
    cnt += 1

# Prepare figure
fig = plt.figure(1)
ax = plt.subplot(111)

for cnt in combinations:
    combination = combinations[cnt]
    files,dates = findFilesWithParameters(metadata, combination)
    n_files = len(files)

    # parameters
    N = int(combination['N'])
    n_steps = int(combination['n_steps'])
    n_save = int(combination['n_save'])
    phi_pack = float(combination['phi_pack'])
    rho=np.sqrt(N/4/phi_pack) # sphere radius
    nu_0 = float(combination['nu_0'])
    eta_n = float(combination['eta_n'])
    J = float(combination['J'])
    dt = float(combination['dt'])
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

    # Plot the result
    line, = ax.plot(t, p_mean, label=r'{0} = {1}'.format(legend_name, combination[variable_compared]))
    color = line.get_color()
    ax.fill_between(t, p_mean-p_std, p_mean+p_std, facecolor=color, alpha=0.3)

ax.legend(loc=4)
plt.xlabel(r'Time ($\tau$)')
plt.ylabel(r'Order parameter')

plt.show()