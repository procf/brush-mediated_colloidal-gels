## Analyze the results of a BD colloid simulation
## NOTE: requires matching Fortran module and solvopt module
## NOTE: this code assumes 1 colloid type (typeid=0)
## NOTE: to select specific analyses, scroll to the bottom
##       of this file and comment out unwanted analyses in the
##       RUN ON A SIMULATION section
"""
## This code performs the following analyses: 
##    * colloid coordination number:
##       - coordination number distribution (Z counts) for each frame
##       - average coordination number (<Z>) for each frame
##    * void size distribution for the final frame OR a selection of frames, 
##      using two methods: 
##       - Torquato’s Pore Size Distribution
##       - Gubbins’s Pore Size Distribution 
##       (requires the solvopt algorithm, included as a separate f90 module)
##    * extracts colloid position data to from GSD to CSV files for all frames
##    * extracts the list of all bonded colloid pairs (i.e network edges) from
##      GSD to CSV files for all frames
##    * calculates primary network analysis metrics for all frames
##       - number of connected components 
##       - average degree (AKA average coordination number in gels)
##       - largest connected component (LCC)
##       - average clustering coefficient
##       - average square clustering coefficient 
"""
## (Rob Campbell)


########################
""" MODULE LIBRARY """
########################
import numpy as np
import pandas as pd
import gsd.hoomd
import math
import networkx as nx
from statistics import mean
from scipy.spatial.distance import squareform
import fortranmod as module
from fortranmod import void_size_calculation
import os
import re
import glob
import sys


##########################
""" INPUT PARAMETERS """
##########################

gel_file = str(sys.argv[1])
tag = str(sys.argv[2])
hb = float(sys.argv[3])
kappa_m = float(sys.argv[4])
d_g = 35 #[nm] = 2*r_g = 2*17.5
R_C_real = 550 #[nm]
colloid_typeid = 0

# filepath to folder where data files will be created
data_outpath = 'data_'+tag

# use brush and depletant to calculate attraction range (cut-off distance)
cut_off = round((d_g+(2*hb))/R_C_real,3)

###########################
""" DEFINE SIM CHECKS """
###########################

# create "data" subfolder if it doesn't exit
if os.path.exists(data_outpath) == False:
  os.mkdir(data_outpath)


#######
## coordination number AKA contact number 
"""
# calculate the coordination number distribution (Z counts) and 
# the average coordination number (<Z>) for each frame
#
# for a given particle, Z is the number of other particles touching it
# (we define "contact" from the "attraction range" set with kappa_m or 
#   brush and depletant parameters size; for a gel the average Z should 
#    plateau over time as the network is formed, and the distribution 
#    is usually centered around an average of Z=6 for central forces)
"""
def coordination_number_py(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename,'r')
  # get the number of frames
  nframes = len(traj)
  # use the last frame to get simulation box size [L_X, L_Y, L_Z]
  Lbox = traj[-1].configuration.box[:3]
  # use the last frame and typeid get the tags and total number colloids
  colloids = np.where(traj[-1].particles.typeid == [colloid_typeid])[0]
  ncolloids = len(colloids)
  # use last frame to get the radius of each colloid
  R_C = 0.5*traj[-1].particles.diameter[colloids]
  # use last frame to get the typeid of each colloid
  typeid = traj[-1].particles.typeid[colloids]

  ## 1. CALCULATE Z DISTRIBUTION

  # gather data for the whole simulation
  allpos_allframe = np.zeros((nframes,ncolloids,3))
  m_xys = np.zeros(nframes)
  for frame in range(nframes):
    # get all particle positions
    allpos_allframe[frame] = traj[frame].particles.position[colloids]
    # get the xy tilt factor (square=0.0, sheared-right=0.45)
    m_xys[frame] = traj[frame].configuration.box[3]

  print('Calculating Z-distribution data for '+str(nframes)+' times...')

  # run the fortran module
  Zs_array = module.coordination_number(nframes,Lbox,ncolloids,R_C,m_xys,allpos_allframe,cut_off)

  # convert array to data frame for easy saving
  allframes_Zs_df = []
  for frame in range(nframes):
    frame_df = pd.DataFrame(Zs_array[frame])
    frame_df.reset_index(inplace=True)
    frame_df = frame_df.rename(columns = {'index':'colloidID'})
    frame_df.insert(loc=1, column='typeid', value=typeid)
    frame_df.insert(loc=0, column='frame', value=frame)
    frame_df = frame_df.rename(columns = {0:'Z'})
    allframes_Zs_df.append(frame_df)

  total_Zsdf = pd.concat(allframes_Zs_df, ignore_index=True)

  total_Zsdf.to_csv(data_outpath+'/Z-counts.csv', index=False)

  print("...Coordination number calculated for all colloids in "+str(nframes)+" frames")


  ## 2. CALCULATE THE AVERAGE Z

  print("Calculating the average coordination number for "+str(nframes)+" frames...")

  # create an array to save the averages
  Zavgs_array = np.zeros(nframes)

  #caluclate the average coordination numbers
  for frame in range(nframes):
    all_sum = sum(Zs_array[frame,:])
    Zavgs_array[frame] = all_sum / ncolloids

  # convert array to data frame for easy saving
  Zavgs_df = pd.DataFrame(Zavgs_array)
  Zavgs_df.insert(loc=0,column='simframe', value=list(range(0,nframes)))
  for frame in range(nframes):
    Zavgs_df = Zavgs_df.rename(columns = {0:'Z_any'})

  Zavgs_df.to_csv(data_outpath+'/Zavg.csv', index=False)

  print("...Average coordination number calculated for "+str(nframes)+" frames")
#######

#######
## void size calculation
"""
# calculates a distribution of the void space (approximated as spheres) in between particle clusters
# 
# There are two methods that are used:
#   - Torquato’s Pore Size Distribution (volume where the center of a particle can fit in between clusters)
#   - Gubbins’s Pore Size Distribution (volume occupied by a whole particlein in between clusters) 
#       (Gubbin's PSD requires solvopt, the Solver For Local Nonlinear Optimization Problems, 
#        incuded as a secondary fortran module) 
# 
# These methods are described in the Section 4.2 and the Appendix of
# Sorichetti, Hugouvieux, and Kob 2020, DOI: 10.1021/acs.macromol.9b02166
#
# This code assumes you are analyzing a porous medium made of uniform particles (like a colloidal gel),
#  and uses particle trajectories in gsd format; It takes particle and probe size as inputs.
#  Then it uses Linked-list method to compute minimum distance of a point in the void space from 
#  nearby porous medium particles. And then uses solvopt non-linear optimization code to compute Gubbin's pore size.
#
# (we use Gubbin's PSD and expect a plot of probability vs void diameter to peak at the size of the most common voids)
"""
# calculate voidsize for the last frame only (or select a different frame)
framechoice_ps = [-1]

def void_size_calc_py(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')
  # get the number of frames
  nframes = len(traj)

  # convert negative framechoice values into specific frames 
  for i in range(len(framechoice_ps)):
    if framechoice_ps[i] < 0:
      framechoice_ps[i] = nframes - abs(framechoice_ps[i])
  # replace the total nframes with the desired nframes to analyze
  nframes = len(framechoice_ps)


  # get the index of all type1 colloids
  colloids = np.where(traj[-1].particles.typeid == [colloid_typeid])[0]
  # calculate the number of type1 colloids
  ncolloids = len(colloids)
  # find the radii of every type1 colloid
  radii = 0.5*traj[-1].particles.diameter[colloids]

  # set the number of random points used to explore void size
  nprobe = 10000 # can test quickly at 1
  # set the size of the cells used in the linked list
  dcell_init = 5.0

  # create empty arrays for holding data for all frames
  box_length = np.zeros((nframes,3))
  rxi = np.zeros((nframes,ncolloids))
  ryi = np.zeros((nframes,ncolloids))
  rzi = np.zeros((nframes,ncolloids))

  # fill arrays with data from each frame
  for i in range(nframes):
    frame = framechoice_ps[i]
    box_length[i,:]=traj[frame].configuration.box[:3]
    rxi[i,:] = traj[frame].particles.position[colloids,0]
    ryi[i,:] = traj[frame].particles.position[colloids,1]
    rzi[i,:] = traj[frame].particles.position[colloids,2]

  # calculate the void_size for all the selected data
  void_size_calculation.void_size_calc(data_outpath,ncolloids,nframes,framechoice_ps,nprobe,radii,dcell_init,rxi,ryi,rzi,box_length)

  print("Pore size distribution calculated for "+str(nprobe)+" probe points in each of "+str(nframes)+" frames")
#######


#######
## account for periodic boundaries in percolation calculation
def pairwise_pbc_dists(positions, box):
    """
    Compute all pairwise distances under PBC using the minimum image convention.
    positions: (N, 3) array
    box: array-like [Lx, Ly, Lz]
    """
    delta = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]  # (N, N, 3)
    delta -= box * np.round(delta / box)  # Apply minimum image convention
    dists = np.linalg.norm(delta, axis=-1)
    return dists
#######


#######
## get cluster span from periodic positions
def unwrap_cluster_positions(positions, box, edges, global_nodes):
    """
    Unwraps positions of a connected cluster so that particles connected by edges
    are placed contiguously in space (removes PBC folding).
    positions: (N,3) positions of cluster nodes
    box: array [Lx, Ly, Lz]
    edges: list of (i,j) pairs defining connectivity within cluster (global IDs)
    global_nodes: list or set of global node IDs corresponding to positions
    """
    # map global → local indices
    global_to_local = {gid: idx for idx, gid in enumerate(global_nodes)}

    # convert edges to local indices
    local_edges = [
        (global_to_local[i], global_to_local[j])
        for (i, j) in edges
        if (i in global_to_local and j in global_to_local)
    ]

    unwrapped = np.copy(positions)
    visited = np.zeros(len(positions), dtype=bool)
    visited[0] = True  # start from first node
    queue = [0]

    while queue:
        i = queue.pop(0)
        for (a, b) in local_edges:
            if a == i or b == i:
                j = b if a == i else a
                if not visited[j]:
                    delta = positions[j] - positions[i]
                    delta -= box * np.round(delta / box)
                    unwrapped[j] = unwrapped[i] + delta
                    visited[j] = True
                    queue.append(j)
    return unwrapped
#######


#######
## posCSV
"""
# create a CSV of particle (i.e. node) position information for each frame
# (tag, x, y, z, typeID, radius)
"""
def posCSV_calc(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')

  # set path and filename
  dir_path = data_outpath+'/frame-pos'
  pos_output = dir_path+'/positions_frame' # + <#>.csv in python loop

  # check for existing CSV data 
  if os.path.exists(dir_path) == False:
    os.mkdir(dir_path)
  if os.path.exists(dir_path) == True:
    # NOTE: this counts ALL CSV files in this directory
    nframes = 0
    # set the pattern for files ending in <number>.csv
    pattern = re.compile(r'positions_frame\d+\.csv$')
    # Iterate directory
    for filename in os.listdir(dir_path):
      # Check if the file has a CSV extension
      if pattern.match(filename):
        nframes += 1

    if nframes == len(traj):
      print('position data CSV files already seem to exist for all frames. Not creating new CSV files.')
      return

  nframes = len(traj)
  colloids = np.where(traj[-1].particles.typeid == [colloid_typeid])[0]
  ncolloids = len(colloids)
  typeid = traj[-1].particles.typeid[colloids]
  radii = 0.5*traj[-1].particles.diameter[colloids]

  for f in range(nframes):
    rpos = traj[f].particles.position[colloids]

    # make a data frame to export to CSV
    df_pos = pd.DataFrame()
    df_pos['tag'] = colloids
    df_pos['x'] = rpos[:,0]
    df_pos['y'] = rpos[:,1]
    df_pos['z'] = rpos[:,2] 
    df_pos['typeID'] = typeid
    df_pos['radius'] = radii

    df_pos.to_csv(pos_output+str(f)+'.csv', index=False)
  print("Position data saved to CSV for "+str(nframes)+" frames")
#######


#######
## edgelistCSV
"""
# create a CSV file of all the bonded particle pairs (i.e. edges of the network) for each frame
# i,j position for each bond/edge
"""

def edgelistCSV_calc(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')
  
  # set path and filename
  edge_dir_path = data_outpath+'/frame-edges'
  edge_output = edge_dir_path+'/edgelist' # + <#>.csv in python loop
  
  # check for existing CSV data 
  if os.path.exists(edge_dir_path) == False:
    os.mkdir(edge_dir_path)
  if os.path.exists(edge_dir_path) == True:
    # NOTE: this counts ALL CSV files in this directory
    nframes = 0
    # set the pattern for files ending in <number>.csv
    pattern = re.compile(r'edgelist\d+\.csv$')
    # Iterate directory
    for filename in os.listdir(edge_dir_path):
      # Check if the file has a CSV extension
      if pattern.match(filename):
        nframes += 1
    
    if nframes == len(traj):
      print('edgelist data CSV files already seem to exist for all frames. Not creating new CSV files.')
      return

  nframes = len(traj)
  colloids = np.where(traj[-1].particles.typeid == [colloid_typeid])[0]
  ncolloids = len(colloids)
  radii = 0.5*traj[-1].particles.diameter[colloids]
  rcut = cut_off
  lbox = traj[-1].configuration.box[:3]

  # create an array of xyz positon of all colloids in all frames    
  allpos = np.zeros((nframes,ncolloids,3))
  for i in range(0,nframes):
    allpos[i,:,:] = traj[i].particles.position[colloids] 
 
  module.edgelist_calc(nframes,ncolloids,radii,allpos,lbox,rcut,edge_output)
  print("Edgelist calculation complete for "+str(nframes)+" frames")
#######


#######
## networkx analysis
"""
# use the networkx package to calculate for all frames:
#   - the number of connected components
#   - average degree (AKA average coordination number)
#   - largest connected component (LCC)
"""

def primary_networkx_calc(filename):

  edge_output = data_outpath+'/frame-edges/edgelist' # + <#>.csv in f90

  # get the number of colloids
  traj = gsd.hoomd.open(filename, 'r')
  nframes = len(traj)
  colloids = np.where(traj[-1].particles.typeid == [colloid_typeid])[0]
  ncolloids = len(colloids)

  # box diagonal
  Lbox = traj[0].configuration.box[:3]
  max_diameter = np.sqrt(Lbox[0]**2 + Lbox[1]**2 + Lbox[2]**2)

  # create empty dataframe for all network data  
  network_frames_df = pd.DataFrame()
  
  # import all data into one dataframe
  frame_dfs = []
  for frame in range(nframes):
    # loop through all frames
    filepath = edge_output+str(frame)+'.csv'
    # import CSV data
    df = pd.read_csv(filepath)
    # rename colums as needed
    df = df.rename(columns={"i": "source", "j": "target"})
    df.insert(loc=0, column='frame', value=frame)
    frame_dfs.append(df)

  alledge_df = pd.concat(frame_dfs, ignore_index=True) 

  # network analysis
  for frame in range(nframes):
    df = alledge_df[alledge_df['frame'] == frame][["source", "target"]]

    # create the network from edge list
    g = nx.from_pandas_edgelist(df)

    # if a node is not in the network, add it
    for particle in range(ncolloids):
      if (  not(   g.has_node(particle)   )  ):
        g.add_node(particle)

    # calculate the total number of edges      
    nedges = nx.number_of_edges(g)

    # calculate the average degree of the network, i.e. avg contact number in gels
    avg_degree = 2 * nedges / ncolloids

    # number of connected components
    n_cc = nx.number_connected_components(g)

    # return the indices of the nodes of the largest connected components
    lcc_nodes = max(nx.connected_components(g), key=len)

    # calculate the size of the largest connected component
    lcc_size = len(lcc_nodes)

    # physical diameter = maximum distance in the lcc
    pos = traj[frame].particles.position[colloids]
    lcc_pos = pos[list(lcc_nodes)]
    if len(lcc_pos) != 0:
      # account for periodic boundaries
      #dists = pairwise_pbc_dists(lcc_pos, Lbox)
      # Only get the upper triangle, since distances are symmetric
      #physical_diameter = np.max(np.triu(dists, k=1))

      # unwrap the LCC before measuring its span
      lcc_edges = list(g.subgraph(lcc_nodes).edges())
      unwrap_pos = unwrap_cluster_positions(lcc_pos, Lbox, lcc_edges, list(lcc_nodes))
      # recentre to remove global drift and align within one periodic image
      unwrap_centered = unwrap_pos - np.min(unwrap_pos, axis=0)
      # compute extent (peak-to-peak) in each direction
      extent = np.ptp(unwrap_centered, axis=0)
      # clip to at most box size in each direction (since periodic)
      extent = np.minimum(extent, Lbox)
      # compute Euclidean span across box
      physical_diameter = np.linalg.norm(extent)

    else:
      physical_diameter = 0.0

    # compile outputs
    data={'frame'                  :[frame],
          'cut_off'                :[cut_off],
          'n_components'           :[n_cc], 
          'lcc_size'               :[lcc_size],
          'ncolloids'              :[ncolloids],
          'lcc_span'               :[physical_diameter],
          'box_span'               :[max_diameter],
          'percolation'            :[physical_diameter/max_diameter],
          'avg_degree'             :[avg_degree],
          }

    res_df = pd.DataFrame(data)
    network_frames_df = pd.concat([network_frames_df, res_df])

  # write all data to one CSV file
  network_frames_df.to_csv(data_outpath+'/networkx-allframes.csv',index = False)
  print("Primary network analysis complete")
#######


###########################
""" RUN ON A SIMULATION """
###########################

if __name__ == '__main__':	
  coordination_number_py(gel_file)
  msd_py(gel_file)
  void_size_calc_py(gel_file)
  posCSV_calc(gel_file)	
  edgelistCSV_calc(gel_file)
  primary_networkx_calc(gel_file)

  print("\nAll analyses complete for:")
  # record params:
  print("\ntag:",tag)
  print("R_C_real:",R_C_real,"[nm]")
  print("hb:",hb,"[nm]")
  print("d_g:",d_g,"[nm]")
  print("kappa_m:",kappa_m)
  print("cut_off:",cut_off)
  print("colloid_typeid:",colloid_typeid)
