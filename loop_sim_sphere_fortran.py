import sim_sphere_fortran
from myutils import readConfigFile
from datetime import datetime
import metadata_lib
import itertools

# ----- PARAMETERS ----- #

parameters = readConfigFile('parameters.ini')

all_N = [100]
all_phi = [1]
all_nu = [1]
all_J = [0.01, 0.05, 0.1, 0.2, 0.4, 0.5, 1, 2]
all_eta = [0.1, 0.5, 1, 2]
all_anisotropy = [1]
all_max_dist = [0]
all_update_nn = [0, 1]

all_comb = [all_N, all_phi, all_nu, all_J, all_eta, all_anisotropy, all_max_dist, all_update_nn]
combinations = list(itertools.product(*all_comb))

# ----- RUN FOR ALL THE PARAMETERS ----- #
repeat = 3
for N, phi_pack, nu_0, J, eta_n, anis, dist, update_nn in combinations:

    # change the parameters
    parameters['N'] = N
    parameters['nu_0'] = nu_0
    parameters['J'] = J
    parameters['eta_n'] = eta_n
    parameters['phi_pack'] = phi_pack
    parameters['fanisotropy'] = anis
    parameters['max_dist'] = dist
    parameters['update_nn'] = update_nn 

    for i in range(0,repeat):
        sim_sphere_fortran.main(parameters)

print('Updating metadata...')
metadata_lib.generateMetadata('./data')
print('Done')