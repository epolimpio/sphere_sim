import numpy as np 
import sim_sphere_fortran as sf

N = 100
r = sf.constrain_to_sphere(2*np.random.rand(3,N)-1,2)
data = np.ones((1,N))
n_bins = 5
print(sf.makePolarData(r, data, n_bins,'sum')[0])