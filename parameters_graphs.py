import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
#import cPickle as pickle # For Python 3 is just pickle
import pickle
from datetime import datetime
import sys

if len(sys.argv)>=2:
    str_date = sys.argv[1]
else:
    str_date = input("Insert date, format Y,m,d,H,M,S ->")

try:
    yy,mm,dd,h,m,s = [int(x) for x in str_date.split(',')[:6]]
except:
    print("ERROR!")
    sys.exit("Invalid date. The format must be Y,m,d,H,M,S separated by comma.")

outfile_analysis=datetime.strftime(datetime(yy,mm,dd,h,m,s),'./data/sphere_data_analysis_%Y-%m-%d_%H-%M-%S.p')

try:
    data=pickle.load(open(outfile_analysis, "rb" ) )
except:
    print("ERROR!")
    sys.exit("Could not read file. Make sure that the file %s exists." % outfile_analysis)

parameters=data[0]
p_angular=data[1]
av_dist_pairs=data[2]

# parameters
N = int(parameters['N'])
n_steps = int(parameters['n_steps'])
n_save = int(parameters['n_save'])
fpacking = float(parameters['fpacking'])
rho=np.sqrt(N/4/fpacking) # sphere radius
F_th_0 = float(parameters['F_th_0'])
dt = float(parameters['dt'])

fig = plt.figure(1)
ax = plt.subplot(111)

ps = [0.0,]*n_steps
for i in range(0,n_steps):
    ps[i] = np.sqrt(p_angular[i].dot(p_angular[i]))/N/rho/F_th_0

t = np.arange(n_steps)*dt
ax.plot(t, ps)
plt.xlabel(r'Time ($\tau$)')
plt.ylabel(r'Order parameter')

fig = plt.figure(2)
ax = plt.subplot(111)

ax.plot(t, av_dist_pairs)
plt.xlabel(r'Time ($\tau$)')
plt.ylabel(r'Average pair distance ($\sigma$)')

plt.show()