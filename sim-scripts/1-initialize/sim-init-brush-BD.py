## Create an initial state init-brush##nm-BD.gsd
## for a Brownian Dynamics simulation with 1 colloid type
## and no overlaps between surface-brush regions
## NOTE: for a variable number of colloids and FIXED BOX SIZE
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
import random # psuedo-random number generator
# Other
import os # miscellaneous operating system interfaces

######### CUSTOM CLASSES
# Function to check overlaps efficiently
def is_valid_position(pos, existing_positions, diameters, curr_diameter, L_X, L_Y, L_Z, tol=0.1):
    if len(existing_positions) == 0:
        return True  # Always accept the first particle

    # Convert list to NumPy array for vector operations
    existing_positions = np.array(existing_positions)

    # Apply periodic boundary conditions
    dr = pos - existing_positions
    dr -= np.round(dr / np.array([L_X, L_Y, L_Z])) * np.array([L_X, L_Y, L_Z])

    # Compute distances
    distances = np.linalg.norm(dr, axis=1)

    # Compute surface-surface distances
    h_ij = distances - 0.5 * (curr_diameter + np.array(diameters[:len(existing_positions)]))

    # Check for overlap
    return np.all(h_ij >= tol)  # Ensure all are at least `tol` apart


#########  SIMULATION INPUTS
# brush size
R_C_nm = 550 # experimental particle size [nm]
h_b_nm = 20  # brush heigh in real units [nm]
h_b = round(h_b_nm/R_C_nm,4) # brush heigh in simulation units

# General parameters
phi = 0.1 # volume fraction
rho = 3   # number density (per unit volume)

# Simulation box size (fixed by L_X)
L_X = 70
L_Y = L_X
L_Z = L_X 
V_total = L_X * L_Y * L_Z # total volume of simulation box (cube)

# Colloid particle details
R_C = 1  # colloid particle radius
V_C = (4./3.) * math.pi * R_C ** 3 # colloid particle volume (1 particle)
m_C = V_C * rho # 1st type colloid particle mass
V_Colloids = phi * V_total # total volume of type 1 colloids
N_C = round(V_Colloids / V_C) # number of 1st type of colloid particles (INT)

# set random seed for repeatable random number generation
seed_value = 42
np.random.seed(seed_value)

######### SIMULATION
## Checks for existing files. If none are found, creates a new 
## random distribution of particles for use as initial simulation state. 

if os.path.exists('init-brush'+str(h_b_nm)+'nm-BD.gsd'):
  print("Initialization file already exists. No new files created.")
else:
  print("New Brownian Dynamics initialization file is being created")
  ## Initialize a snapshot of the system
  snapshot = gsd.hoomd.Frame()
  snapshot.configuration.box = [L_X, L_Y, L_Z, 0, 0, 0] # create the sim box
  snapshot.particles.N=N_C # add all particles to the snapshot
  # set the particle types
  snapshot.particles.types = ['B']
  # assign particles to each type
  typeid = []
  typeid.extend([0]*N_C)
  snapshot.particles.typeid = typeid
  # set a mass for each particle type
  mass = []
  mass.extend([m_C]*N_C)
  snapshot.particles.mass = mass
  # set a diameter for each particle type
  diameter = []
  brush_diameter = []
  diameter.extend([2.0*R_C]*N_C)
  brush_diameter = list(np.array(diameter) + (2*h_b))
  snapshot.particles.diameter = diameter
  # randomly distribute all the particles in 3D space
  # ensuring the brush regions do not overlap
  pos_arr = np.zeros((N_C,3))
  print("\n  Distributing brush coated particles (no brush-brush overlap)")
  sorted_indices = np.argsort(-np.array(brush_diameter))
  placed_particles = 0
  for curr_index in sorted_indices:
    overlap = True
    while overlap:
      new_pos = np.random.uniform(-0.5*np.array([L_X, L_Y, L_Z]), 0.5*np.array([L_X, L_Y, L_Z]))
      if is_valid_position(new_pos, pos_arr[:placed_particles], brush_diameter, brush_diameter[curr_index], L_X, L_Y, L_Z):
        pos_arr[placed_particles] = new_pos
        placed_particles += 1
        overlap = False
  print(f"  ...{placed_particles} particles placed.\n")
  snapshot.particles.position = pos_arr

  # save the snapshot of the initialized system
  with gsd.hoomd.open(name='init-brush'+str(h_b_nm)+'nm-BD.gsd', mode='w') as f:
    f.append(snapshot)

  print("New Brownian dynamics initialization file (init-BD.gsd) created.\n")
  print("Seed for Random Number Generator: "+str(seed_value)) 
  print("Simulation volume: L_X = " + str(L_X) + ", L_Y = " + str(L_Y) 
    + ", L_Z = " + str(L_Z))
  print("Volume fraction: " + str(phi) + "\n")
  print("Total number of colloid particles: " + str(N_C))
  print("Brush length: " + str(h_b_nm) + "nm")
