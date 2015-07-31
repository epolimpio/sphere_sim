import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
#import cPickle as pickle # For Python 3 is just pickle
import pickle

infile='neighbors.p'

data=pickle.load(open( infile, "rb" ) )

number_of_neighbors=np.array(data[0])
number_of_new_neighbors = np.array(data[1])
histogram = np.array(data[2])

print(number_of_new_neighbors/number_of_neighbors)

N = 100
L=0.4
R=0.4
rho=np.sqrt(N*L**2/(4*np.pi))
dr = 0.04
r_vec = dr*np.arange(np.floor(rho/dr)+1)


ax = plt.subplot(111)
ax.plot(histogram[1:]/r_vec[1:]**2)
plt.show()
