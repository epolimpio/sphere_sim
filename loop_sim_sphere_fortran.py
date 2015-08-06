import sim_sphere_fortran
from myutils import readConfigFile
from datetime import datetime

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

all_N = [50, 100, 200]
all_phi = [1]
all_nu = [0.5, 1, 2]
all_J = [0.5, 1, 2]
all_eta = [0.5, 1, 2]

combinations = []
for N in all_N:
    for phi in all_phi:
        for nu in all_nu:
            for J in all_J:
                for eta in all_eta:
                    combinations.append((N,phi,nu,J,eta))

# ----- RUN FOR ALL THE PARAMETERS ----- #
repeat = 3
for N, phi_pack, nu_0, J, eta_n in combinations:

    # change the parameters
    parameters['N'] = N
    parameters['nu_0'] = nu_0
    parameters['J'] = J
    parameters['eta_n'] = eta_n
    parameters['phi_pack'] = phi_pack

    for i in range(0,repeat):
        sim_sphere_fortran.main(parameters)