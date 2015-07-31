import numpy as np 
from fortran_lib import simforces
from math import pi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri

# theta = pi*np.random.rand(N) - pi/2
# phi = 2*pi*np.random.rand(N) - pi

# x = np.zeros(N)
# y = np.zeros(N)
# z = np.zeros(N)

# stripack.trans(theta, phi, x, y, z)

# x = 2*np.random.rand(N) - 1
# y = 2*np.random.rand(N) - 1
# z = 2*np.random.rand(N) - 1

# r = np.sqrt(x**2+y**2+z**2)

# x = x/r
# y = y/r
# z = z/r

# dist = np.zeros(N)
# iwk = np.zeros(2*N)
# list_,lptr,lend,lnew,near,next,dist,ier = stripack.trmesh(x,y,z,iwk[0:N],iwk[N:],dist)
# if ier == 0:
# 	ltri = np.zeros((3,2*N-4))
# 	nt, ltri, ier = stripack.trlist2(list_,lptr,lend,ltri)

# 	if ier == 0:
# 		print(ltri)
# 		fig = plt.figure()
# 		ax = fig.add_subplot(1, 1, 1, projection='3d')
# 		ax.plot_trisurf(x, y, z, triangles=ltri.T-1, cmap=plt.cm.Spectral, alpha=0.5)
# 		plt.show()

print(simforces.calc_force_elastic.__doc__)