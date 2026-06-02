## MPI simulation for Brownian Dynamics gelation of an
## attractive colloid suspension with brush-mediated
## angular constraints on particle motion
## see: https://arxiv.org/abs/2603.13596 
## 
## NOTE: requires matching init-brush##nm-BD.gsd file
## (Rob Campbell)


######### MODULE LIBRARY 
# Use HOOMD-blue
import hoomd
import hoomd.md # molecular dynamics
# Use GSD files
import gsd # provides Python API; MUST import sub packages explicitly (as below)
import gsd.hoomd # read and write HOOMD schema GSD files
# Maths
import numpy as np
import math
import random # pseudo-random number generator
# Other
import os # miscellaneous operating system interfaces
import re # for regex identification of params from path
import sys # pass variables via system from jobscript

#pip install mpi4py
from mpi4py import MPI # for tracking MPI processes

# track MPI rank and comm.size (1 if serial)
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if comm.size > 1:
  print("MPI rank: "+str(rank))


######### SET FILETAG FROM jobscript
tag = str(sys.argv[1])
R_C_real = 550 # real colloid radius [nm]

######### SIMULATION INPUTS
# General parameters
phi =  0.1 # volume fraction
rho = 3.0  # number density (per unit volume)
kT = 0.1   # system temperature
m_K = float(sys.argv[9])

# MorseBrush attraction potential parameters
D0 = float(sys.argv[2]) * kT # attraction strength 
kappa_m = float(sys.argv[3]) # range of attraction ; distance in BD units is approx 3/kappa
r0 = round(float(sys.argv[7])/R_C_real,4) # h_ij position of the depletion minimum [sim units]
hb = round(float(sys.argv[8])/R_C_real,4) # polymer brush length [sim units]
#d_p = 3/kappa_m
d_p = 2 * hb # diameter of the polymer brush [sim units]

# Electrostatic repulsion parameters
kappa_e = float(sys.argv[4])  # inverse Debye length
hb_nm = float(sys.argv[8])    # [nm]
Z = float(sys.argv[5]) * kT   # Z_{H} ~ 2*Z_{L} 
L_elec = float(sys.argv[6])   # [nm]
e0 = round(L_elec/R_C_real,4) 


N_time_steps = 90000000 # number of time steps
dt_Integration = 1e-4   # timestep size
period = 900000         # recording interval [timesteps]


# Colloid particle details
R_C = 1 # colloid particle radius

# Brownian parameters
eta0 = 1.0 # viscosity of the fluid (tunable parameter, not direct viscosity)
gamma = 6.0*np.pi*eta0*R_C # BD stock friction coefficient

# Particle interaction parameters
r_c = 1.0 # cut-off radius parameter, r_c>=3/kappa (r_cut = # * r_c) 
if r_c < (3/kappa):
  if rank == 0 or comm.size == 1:
    print('WARNING: r_c is less than range of attraction. Increase r_c')

seed_value = 42

######### SIMULATION
# bring an equilibrium state towards a quasi-steady state gel

if rank == 0 or comm.size == 1:
  print('Brownian Dynamics equilibrium state is being brought to quasi-steady state gelation')

## Create a CPU simulation
device = hoomd.device.CPU()
# set seed to a fixed value for reproducible simulations
sim = hoomd.Simulation(device=device, seed=seed_value) 

# start the simulation from Equilibrium (and don't set timestep to 0)
sim.create_state_from_gsd(filename=f'/projects/props/Rob/colloids/bimodal/r1-0/BD/L70/mono/seed42/phi10/1-initialize/init-brush{hb_nm}nm-BD.gsd')


# assign particle types to groups
# (in case we want to integrate over subpopulations only,
# but would require other mods to source code)
groupA = hoomd.filter.Type(['B'])
all_ = hoomd.filter.Type(['B'])

# don't want to thermalize the system after equilibrium because
# we are tracking velocity during gelation
#sim.state.thermalize_particle_momenta(filter=all_, kT=kT)

# create neighboring list
nl = hoomd.md.nlist.Tree(buffer=0.05);

# define MorseBrush and Electrostatic interactions
morsebrushe0 = hoomd.md.pair.MorseBrushE0(nlist=nl, default_r_cut=1.0 * r_c, K = m_K, distance = d_p)

# colloid-colloid: hard particles (no deformation/overlap)
morsebrushe0.params[('B','B')] = dict(D0=D0, alpha=kappa_m, r0=r0, soft_shift=1.0, scaled_D0=False, kappa=kappa_e, Z=Z, e0=e0, a1=R_C, a2=R_C)
morsebrushe0.r_cut[('B','B')] = r_c+(R_C+R_C) # used to assemble nl


# choose integration method for end of each timestep | BROWNIAN (overdamped) 
brownian = hoomd.md.methods.Brownian(filter=all_, kT=kT, default_gamma=gamma)
integrator=hoomd.md.Integrator(dt=dt_Integration, forces=[morsebrushe0], methods=[brownian])

sim.operations.integrator = integrator

# set the simulation to log certain values
logger = hoomd.logging.Logger()
thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=all_)
sim.operations.computes.append(thermodynamic_properties)
logger.add(thermodynamic_properties,quantities=['kinetic_temperature',
  'pressure_tensor','virial_ind_tensor','potential_energy'])
logger.add(sim, quantities=['tps'])

# set output file
filename = "Gelation-" + tag + "-BD.gsd"

# Check if the file already exists; if not, proceed with writing
if not os.path.exists(filename):
    gsd_writer = hoomd.write.GSD(trigger=period, filename=filename, 
                                 filter=all_, mode='wb', dynamic=['property', 'momentum', 'attribute'])
    gsd_writer.write_diameter = True
    sim.operations.writers.append(gsd_writer)
    gsd_writer.logger = logger

    if rank == 0 or comm.size == 1:
        # Record simulation parameters
        print("\nR_C:", R_C_real, 'nm (', R_C,')')
        print("phi:", phi)
        print("seed_value:", seed_value,'\n')
        print("Depletion/Morse:")
        print(" - D0/kT:", round(D0 / kT))
        print(" - kappa_m:", kappa_m)
        print(" - r0:", float(sys.argv[7]), 'nm (', r0, ')\n')
        print("Electrostatics:")
        print(" - inverse Debye length:", kappa_e)
        print(" - charge constant Z/(kT*R_C): ", Z/kT)
        print(" - thickness of the particle 'charge layer':", L_elec, 'nm (', e0, ')\n')
        print("Brush Effects:")
        print(" - brush size (hb):", hb_nm, 'nm (', hb, ')')
        print(" - brush rigidity (K):", m_K, '\n')

    # Run simulation
    n_periods = int(np.ceil(N_time_steps / period))
    for i in range(n_periods):
        if i == 0:
            # Write the initial state (e.g., the last frame of Equilibrium) to the gel file
            sim.run(period, write_at_start=True)
        else:
            sim.run(period)
        gsd_writer.flush()

    if rank == 0 or comm.size == 1:
        print(f'New MorseBrush gelation state ({filename}) created.')

else:
    print(f"File {filename} already exists. Skipping simulation.")
