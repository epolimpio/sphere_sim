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

# Parameters
files_to_compare = [
'2015,08,05,21,23,04',
'2015,08,05,22,33,24',
'2015,08,05,23,39,54',
]

variable_compared = 'nu_0'
legend_name = r'$k^{-1}$'

# Prepare figure
fig = plt.figure(1)
ax = plt.subplot(111)

rotation = []
for str_date in files_to_compare:
    yy,mm,dd,h,m,s = [int(x) for x in str_date.split(',')[:6]]
    
    # Open analysis file
    analysis_filename = './data/sphere_data_analysis_%Y-%m-%d_%H-%M-%S.p'
    dtime_data = datetime(yy,mm,dd,h,m,s)
    data = openPickleFile(analysis_filename, dtime_data)

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
    eta_n = float(parameters['eta_n'])
    J = float(parameters['J'])
    dt = float(parameters['dt'])

    print('N: {0}, phi: {1}, J: {2}, nu: {3}, eta: {4}'.format(N, phi_pack, J, nu_0, eta_n))

    ps = [0.0,]*n_steps
    for i in range(0,n_steps):
        ps[i] = np.sqrt(p_angular[i].dot(p_angular[i]))/N/rho/nu_0

    t = np.arange(n_steps)*dt
    
    line, = ax.plot(t, ps, label=r'{0} = {1}'.format(legend_name, parameters[variable_compared]))
    color = line.get_color()

    # If rotated, plot rotated data
    if rotated:
        rotated_filename = str(parameters['outfile_postrotation'])
        data = openPickleFile(rotated_filename, dtime_data)
        rotation.append((data[1]*dt, color))



if rotation:        
    y1, y2 = ax.get_ylim()
    print('Rotation times')
    for rot_time, color in rotation:
        print(rot_time)
        ax.plot([rot_time, rot_time], [y1, y2], '--', color=color) 

ax.legend(loc=4)
plt.xlabel(r'Time ($\tau$)')
plt.ylabel(r'Order parameter')

plt.show()