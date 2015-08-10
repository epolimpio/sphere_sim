import sim_sphere_fortran
from myutils import readConfigFile
from datetime import datetime
import metadata_lib
import itertools

# ----- PARAMETERS ----- #

parameters = readConfigFile('parameters.ini')

# simulation parameters
n_steps = int(parameters['n_steps'])
n_save = int(parameters['n_save'])
ftype = int(parameters['ftype'])

# timestep
dt = float(parameters['dt'])

# outfile for data
tnow = datetime.now()
outfile_video = datetime.strftime(tnow,str(parameters['outfile_video']))
outfile_analysis = datetime.strftime(tnow,str(parameters['outfile_analysis']))
outfile_postrotation = datetime.strftime(tnow,str(parameters['outfile_postrotation']))

# ----- PARAMETERS TO LOOP ----- #
# We write the parameters as tuples (N, phi_pack, nu_0, J, eta_n)

all_N = [100]
all_phi = [1]
all_nu = [0.01, 0.05, 0.1, 0.2, 0.4, 0.5]
all_J = [0.01, 0.05, 0.1, 0.2, 0.4, 0.5]
all_eta = [0.5, 1, 2]
all_anisotropy = [0.5, 1, 2]
all_max_dist = [0, 2.5, 3]

all_comb = [all_N, all_phi, all_nu, all_J, all_eta, all_anisotropy, all_max_dist]
combinations = list(itertools.product(*all_comb))

# ----- RUN FOR ALL THE PARAMETERS ----- #
repeat = 3
for N, phi_pack, nu_0, J, eta_n, anis, dist in combinations:

    # change the parameters
    parameters['N'] = N
    parameters['nu_0'] = nu_0
    parameters['J'] = J
    parameters['eta_n'] = eta_n
    parameters['phi_pack'] = phi_pack
    parameters['fanisotropy'] = anis
    parameters['max_dist'] = dist 

    for i in range(0,repeat):
        sim_sphere_fortran.main(parameters)

print('Updating metadata...')
metadata_lib.generateMetadata('./data')
print('Done')