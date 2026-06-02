## Build particle networks for GMM clustering 
## NOTE: requires matching Fortran module
## NOTE: this code assumes 1 colloid type (typeid=0)
##
"""
## Load data and build a network
##
## How to use:
##     python sim-build-particle-network.py brush_density file_type
##
## - Load data from gel_file (GSD) and build a network, G
## - Collect angle distribution for G
##
## OUTPUT:
##    - node_df_all, node_df_lcc, 
##        # keys : 'Tag#', 'x', 'y', 'z', 'Radius', 'TypeID'  
##        saved to data_outpath+'/node_df_all.csv' 
##        saved to data_outpath+'/node_df_lcc.csv' 
##    - edge_df_all, edge_df_lcc
##        # keys : 'source', 'target', 'edge_type'  
##        saved to data_outpath+'/edgelist_all.csv'
##        saved to data_outpath+'/edgelist_lcc.csv'
##    - G, g
##        # keys : 'pos' [x,y,z], 'radius', 'type_id', 'edge_type' 
##        saved to data_outpath+'/particle-network_G.pkl'
##        saved to data_outpath+'/lcc-particle-network_g.pkl'
##    - lcc_df
##        # keys : 'cut_off', 'n_components', 'lcc_size', 'ncolloids', 'avg_degree', 
##        saved to data_outpath+'/lcc.csv'
##    - angle_df, lcc_angle_df 
##        # keys : 'angle_id', 'edge_1', 'type_1', 'edge_2', 'type_2', 'angle_degree' 
##        saved to data_outpath+'/angle_dist_all.csv'
##        saved to data_outpath+'/angle_dist_lcc.csv'
"""
## (Rob Campbell)

########################
""" MODULE LIBRARY """
########################
# load data and build a network
import numpy as np
import gsd.hoomd
import fortranmod as module 
import pandas as pd
import networkx as nx
import os
import re
import glob
import sys
import pickle
# for plotting
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
# angle gradient matrix
import scipy.linalg as la

##########################
""" INPUT PARAMETERS """
##########################

brush_density = str(sys.argv[1]) # 'H' high 'L' low

if brush_density == 'H':
  gel_file = '/sim-scripts/potential-morsebrush/K_1xUC/hb_20nm-H/6mgmL-10mM-K10/Gelation-6mgmL-10mM-K10_H-BD.gsd'
  tag ='6mgmL-10mM-H_sim' 
  hb = 20   
  kappa_m = 75 
  m_K = 10
  rmin = 12 # [nm]
  U0 = 9.7 

if brush_density == 'L':
  gel_file = '/sim-scripts/potential-morsebrush/K_1xUC/hb_20nm-L/6mgmL-10mM-K10/Gelation-6mgmL-10mM-K10_L-BD.gsd'
  tag = '6mgmL-10mM-L_sim'
  hb = 20
  kappa_m = 75 
  m_K = 20
  rmin = 4 # [nm]
  U0 = 20.4 # * kT

# real parameters
d_g = 35 #[nm] = 2*r_g = 2*17.5
R_C_real = 550 #[nm]

# simulation parameters
R_C = 1
PBC = True
cut_off = round((d_g+(2*hb))/R_C_real,3)
colloid_typeid = 0
kT = 0.1
V_colloid = (4/3)*np.pi*R_C**3

# filepath to folder where data files will be created
data_outpath = 'data_'+tag

# create "data" subfolder if it doesn't exit
if os.path.exists(data_outpath) == False:
  os.mkdir(data_outpath)


###############################
""" BUILD NETWORK FROM DATA """
###############################



# function to add edge lengths to the graph:
def add_edge_lengths(G, box=None):
    """
    Compute physical edge lengths for all edges in a graph G, using node 'pos'.
    Optionally apply periodic boundaries.

    Parameters
    ----------
    G : networkx.Graph
        Graph where each node has attribute 'pos' = (x, y, z)
    box : array-like or None
        If provided, should be shape (3,) = [Lx, Ly, Lz].
        Uses minimum-image convention for periodic boundary conditions.

    Returns
    -------
    None (modifies G in place, adds edge attribute 'length')
    """

    use_pbc = box is not None
    if use_pbc:
        box = np.asarray(box, dtype=float)

    for u, v, data in G.edges(data=True):
        # positions from node attributes
        r1 = np.asarray(G.nodes[u]['pos'], dtype=float)
        r2 = np.asarray(G.nodes[v]['pos'], dtype=float)

        delta = r2 - r1

        if use_pbc:
            # minimum-image
            delta -= box * np.round(delta / box)

        length = np.linalg.norm(delta)

        # write edge attribute
        data['length'] = length

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
    for fn in os.listdir(dir_path):
      # Check if the file has a CSV extension
      if pattern.match(fn):
        nframes += 1

    if nframes == len(traj):
      print(' - position data CSV files already seem to exist for all frames. Not creating new CSV files.')
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
  print(" - Position data saved to CSV for "+str(nframes)+" frames")

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
    for fn in os.listdir(edge_dir_path):
      # Check if the file has a CSV extension
      if pattern.match(fn):
        nframes += 1
  
    if nframes == len(traj):
      print(' - edgelist data CSV files already seem to exist for all frames. Not creating new CSV files.')
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
  print(" - Edgelist calculation complete for "+str(nframes)+" frames")

# create the edgelist
posCSV_calc(gel_file)	
edgelistCSV_calc(gel_file)

# extract additional data from the GSD file 
traj = gsd.hoomd.open(gel_file, 'r')
nframes = len(traj)
colloids = np.where(traj[-1].particles.typeid == [colloid_typeid])[0]
ncolloids = len(colloids)
pos = traj[-1].particles.position[colloids]
radius = 0.5*traj[-1].particles.diameter[colloids]
typeid = traj[-1].particles.typeid[colloids]
Lbox = traj[0].configuration.box[:3]

system_volume = Lbox[0] * Lbox[1] * Lbox[2]
data = {
  "Tag#": colloids,
  "x": pos[:,0],
  "y": pos[:,1],
  "z": pos[:,2],
  "Radius": radius,
  "TypeID": typeid,
}
df = pd.DataFrame(data)

# create empty dataframe for all network data  
network_frames_df = pd.DataFrame()
  
# import all data into one dataframe
edge_output = data_outpath+'/frame-edges/edgelist' # + <#>.csv in f90
frame_dfs = []
for frame in range(nframes):
  # loop through all frames
  filepath = edge_output+str(frame)+'.csv'
  # import CSV data
  edge_df = pd.read_csv(filepath)
  # rename colums as needed
  edge_df = edge_df.rename(columns={"i": "source", "j": "target"})
  edge_df.insert(loc=0, column='frame', value=frame)
  frame_dfs.append(edge_df)

alledge_df = pd.concat(frame_dfs, ignore_index=True) 

# only look at the last frame 
edge_df = alledge_df[alledge_df['frame'] == (nframes-1)][["source", "target"]]

#######
## CREATE A NEW GRAPH
G = nx.Graph()

# track the type of particle-particle interaction in each edge
edge_id = []
for idx, edge_data in edge_df.iterrows():
    source_tag = edge_data['source']
    source_type_id = df.loc[df['Tag#'] == source_tag, 'TypeID'].values[0]

    target_tag = edge_data['target']
    target_type_id = df.loc[df['Tag#'] == target_tag, 'TypeID'].values[0]

    # Create the edge_type by combining source and target type_ids
    edge_type = f"{source_type_id}-{target_type_id}"
    edge_id.append(edge_type)

# Add the new 'edge_type' column to edge_df
edge_df['edge_type'] = edge_id

# Add nodes to the graph from edge_df (Tag# as the node ID, and positions as attributes)
for _, row in df.iterrows():
    G.add_node(row['Tag#'], pos=(row['x'], row['y'], row['z']), radius=row['Radius'], type_id=row['TypeID'])

# Add edges to the graph from edge_df (target and source define the edge connections, edge_type as attribute)
for _, row in edge_df.iterrows():
    G.add_edge(row['target'], row['source'], edge_type=row['edge_type'])

add_edge_lengths(G, box=Lbox) # with PBC
# Extract all edge lengths from the Graph
rows = []
for u, v, data in G.edges(data=True):
    rows.append({
        "source": u,
        "target": v,
        "edge_length": data.get("length", np.nan)
    })
edge_len_df = pd.DataFrame(rows)
# Merge the lengths back into edge_df
# Note: Your edge_df uses ['target','source'], so match on both
edge_df = edge_df.merge(edge_len_df,
                       on=["source", "target"],
                       how="left")


# Access the node positions
pos = nx.get_node_attributes(G, 'pos')

print(' - Full network:')

# save as csv files
df.to_csv(data_outpath+'/node_df_all.csv',index=False) 
print('   - GSD data saved as CSV: "node_df_all.csv"')

edgelist_file = data_outpath+'/edgelist_all.csv'
edge_df.to_csv(edgelist_file,index=False) 
print('   - Edgelist calculated, saved to "edgelist_all.csv"')

# SAVE GRAPH AS PKL FILE
particle_network_filename = "particle-network_G.pkl"
with open(data_outpath+'/'+particle_network_filename, 'wb') as fr:
  pickle.dump(G, fr)
print('   - Graph saved to "particle-network_G.pkl" (PKL file)')

phi = (ncolloids*V_colloid) / system_volume

# calculate the total number of edges      
nedges = nx.number_of_edges(G)

# calculate the average degree of the network, i.e. avg contact number in gels
avg_degree = 2 * nedges / ncolloids

# number of connected components
n_cc = nx.number_connected_components(G)

# return the indices of the nodes of the largest connected components
lcc_nodes = max(nx.connected_components(G), key=len)

# calculate the size of the largest connected component
lcc_size = len(lcc_nodes)

# compile outputs
data={
      'cut_off'                :[cut_off],
      'n_components'           :[n_cc], 
      'lcc_size'               :[lcc_size],
      'ncolloids'              :[ncolloids],
      'avg_degree'             :[avg_degree],
      'phi'                    :[phi],
      'system_volume'          :[system_volume],
      'L_X'                    :[Lbox[0]],
      'L_Y'                    :[Lbox[1]],
      'L_Z'                    :[Lbox[2]],
      }

net_df = pd.DataFrame(data)


###################################
""" FUNCTION TO UNWRAP POSITIONS"""
###################################

def unwrap_cluster_positions(positions_df, box, edges, global_nodes, start_nodes=None):
    """
    Unwraps positions for the given node set (global_nodes) *and* the provided edges.
    This function unwraps each connected component inside the provided node set separately.
    Returns ndarray shape (n_nodes, 3) in the same order as global_nodes.
    """
    # map global Tag# → local index
    global_to_local = {gid: idx for idx, gid in enumerate(global_nodes)}

    # local edges (only those that connect two nodes in this cluster)
    local_edges = [
        (global_to_local[i], global_to_local[j])
        for (i, j) in edges
        if (i in global_to_local and j in global_to_local)
    ]

    # reorder positions according to global_nodes -> shape (n,3)
    pos_df_indexed = positions_df.set_index("Tag#")
    positions = pos_df_indexed.loc[global_nodes][["x","y","z"]].values.astype(float)

    n = len(positions)
    unwrapped = positions.copy()
    visited = np.zeros(n, dtype=bool)

    # build adjacency list for local indices
    adj = [[] for _ in range(n)]
    for a, b in local_edges:
        adj[a].append(b)
        adj[b].append(a)

    # find connected components in the local graph (indices)
    comps = []
    for idx in range(n):
        if not visited[idx]:
            # BFS/DFS to get this component
            comp_nodes = []
            stack = [idx]
            visited[idx] = True
            while stack:
                u = stack.pop()
                comp_nodes.append(u)
                for v in adj[u]:
                    if not visited[v]:
                        visited[v] = True
                        stack.append(v)
            comps.append(comp_nodes)

    # For each component, perform unwrapping starting from an arbitrary node in that component
    for comp in comps:
        # pick a root
        root = comp[0]
        comp_visited = {root}
        queue = [root]
        while queue:
            i = queue.pop(0)
            for j in adj[i]:
                if j not in comp_visited:
                    # minimum-image delta using box
                    delta = positions[j] - positions[i]
                    delta -= box * np.round(delta / box)
                    unwrapped[j] = unwrapped[i] + delta
                    comp_visited.add(j)
                    queue.append(j)
        # At this point, the component nodes are internally unwrapped relative to each other

    return unwrapped

###############################
""" GET LCC DATA SEPARATELY """
# these systems are fully percolated, so full network and LCC should be the same
# however, for data that is not fully-connected we expect the LCC to be the mechanically
# relevant portion of the network; this can be compared to rattler-removal based on
# the isostaticity criterion. In our systems the results from the LCC, and rattler-removal
# are the same. Here we use the LCC for internal consistency
###############################

# (a) calculate the LCC (largest connected component)
#     this is the part of the network that contributes to ELASTICITY
g = G.subgraph(lcc_nodes).copy() # Create a subgraph containing only the lcc

# filter data for lcc only (used in plotting)
node_df_lcc = df[df['Tag#'].isin(lcc_nodes)].copy()

edge_df_lcc = edge_df[(edge_df['target'].isin(lcc_nodes)) & (edge_df['source'].isin(lcc_nodes))].copy()
ncolloids_lcc = len(node_df_lcc)
print(' - LCC contains',ncolloids_lcc,'colloids ('+str(round(ncolloids_lcc/ncolloids,4)*100)+'% of all colloids)')

# save lcc data for easier access later
lcc_node_file = data_outpath+'/node_df_lcc.csv'
node_df_lcc.to_csv(lcc_node_file,index=False)
lcc_edgelist_file = data_outpath+'/edgelist_lcc.csv'
edge_df_lcc.to_csv(lcc_edgelist_file,index=False)
print('   - LCC data saved to "node_df_lcc.csv", "edgelist_lcc.csv"')

# SAVE GRAPH AS PKL FILE
lcc_particle_network_filename = "lcc-particle-network_g.pkl"
with open(data_outpath+'/'+lcc_particle_network_filename, 'wb') as fr:
  pickle.dump(g, fr)
print('   - LCC graph saved to "lcc-particle-network_g.pkl" (PKL file)')

phi_lcc = (ncolloids_lcc*V_colloid) / system_volume

# calculate the total number of edges      
nedges_lcc = nx.number_of_edges(g)

# calculate the average degree of the network, i.e. avg contact number in gels
avg_degree_lcc = 2 * nedges_lcc / ncolloids_lcc

# compile outputs
data={
      'cut_off'                :[cut_off],
      'ncolloids'              :[ncolloids_lcc],
      'avg_degree'             :[avg_degree_lcc],
      'phi'                    :[phi_lcc],
      'system_volume'          :[system_volume],
      'L_X'                    :[Lbox[0]],
      'L_Y'                    :[Lbox[1]],
      'L_Z'                    :[Lbox[2]],
      }

lcc_df = pd.DataFrame(data)



##########################
""" COLLECT ANGLE DATA """
##########################


# Function to compute the angle between two vectors in 3D
def angle_between_vectors(v1, v2):
    # Normalize vectors
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    # Dot product and calculate the angle
    dot_product = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
    angle = np.arccos(dot_product)
    return np.degrees(angle)

def pbc_vec(a, b, box):
    dv = b - a
    return dv - box * np.round(dv / box)


## FULL NETWORK

# Prepare a list to store angle data
angle_data = []

# Angle ID counter
angle_id_counter = 0

# ensure node labels are integers
G = nx.relabel_nodes(G, lambda x: int(x))

# Retrieve node positions and edge types
pos = nx.get_node_attributes(G, 'pos')
edge_types = nx.get_edge_attributes(G, 'edge_type')

# Calculate the angle between edges at each node
angle_dict = {}
for node in G.nodes():
    neighbors = list(G.neighbors(node))
    if len(neighbors) < 2:
        continue  # No angle if less than two edges meet at the node
    
    for i in range(len(neighbors)):
        for j in range(i + 1, len(neighbors)):
            # Get position vectors for two neighbors
            if PBC == True:
              vec1 = pbc_vec(np.array(pos[node]), np.array(pos[neighbors[i]]), Lbox) # PBC
              vec2 = pbc_vec(np.array(pos[node]), np.array(pos[neighbors[j]]), Lbox) # PBC
            else:
              vec1 = np.array(pos[neighbors[i]]) - np.array(pos[node])
              vec2 = np.array(pos[neighbors[j]]) - np.array(pos[node])
            
            # Calculate the angle between the two vectors
            angle = angle_between_vectors(vec1, vec2)

            # Store the angle as a node attribute (angles between edges meeting at this node)
            G.nodes[node][f'angle_{neighbors[i]}_{neighbors[j]}'] = angle
            
            # Store the angle between the two neighbors at the node in a dictionary
            angle_dict[(node, neighbors[i], neighbors[j])] = angle
            
            # Retrieve edge types for the two edges
            edge_1 = (node, neighbors[i]) if (node, neighbors[i]) in edge_types else (neighbors[i], node)
            edge_2 = (node, neighbors[j]) if (node, neighbors[j]) in edge_types else (neighbors[j], node)
            
            type_1 = edge_types[edge_1]
            type_2 = edge_types[edge_2]
            
            # Increment angle_id
            angle_id_counter += 1
            
            # Store the information in a list of tuples
            angle_data.append((angle_id_counter, edge_1, type_1, edge_2, type_2, angle))

# Create a DataFrame from the angle data
angle_df = pd.DataFrame(angle_data, columns=["angle_id", "edge_1", "type_1", "edge_2", "type_2", "angle_degree"])

## LCC angles

# ensure node labels are integers
g = nx.relabel_nodes(g, lambda x: int(x))

# Prepare a list to store angle data
lcc_angle_data = []

# Retrieve node positions and edge types
#lcc_pos = nx.get_node_attributes(g, 'pos')
unwrapped_pos_arr = unwrap_cluster_positions(node_df_lcc[['Tag#','x','y','z']], Lbox, edges=list(g.edges()), global_nodes=list(g.nodes()))
nodes = list(g.nodes())
# Convert array → dict mapping node → (x,y,z)
unwrapped_pos_lcc = {
    nodes[i]: unwrapped_pos_arr[i]
    for i in range(len(nodes))
}
for n in g.nodes():
    g.nodes[n]['pos'] = unwrapped_pos_lcc[n]
lcc_pos = nx.get_node_attributes(g, 'pos')
unwrapped_df = pd.DataFrame(
    unwrapped_pos_arr,
    index=nodes,
    columns=['x_unwrapped', 'y_unwrapped', 'z_unwrapped']
)
unwrapped_df.index.name = 'Tag#'
node_df_lcc = (
    node_df_lcc
    .set_index('Tag#')
    .join(unwrapped_df, how='left')
    .reset_index()
)

lcc_edge_types = nx.get_edge_attributes(g, 'edge_type')

# Angle ID counter
lcc_angle_id_counter = 0

# Calculate the angle between edges at each node
lcc_angle_dict = {}
for node in g.nodes():
    neighbors = list(g.neighbors(node))
    if len(neighbors) < 2:
        continue  # No angle if less than two edges meet at the node

    for i in range(len(neighbors)):
        for j in range(i + 1, len(neighbors)):
            # Get position vectors for two neighbors
            if PBC == True:
              vec1 = pbc_vec(np.array(lcc_pos[node]), np.array(lcc_pos[neighbors[i]]), Lbox) # PBC
              vec2 = pbc_vec(np.array(lcc_pos[node]), np.array(lcc_pos[neighbors[j]]), Lbox) # PBC
            else:
              vec1 = np.array(lcc_pos[neighbors[i]]) - np.array(lcc_pos[node])
              vec2 = np.array(lcc_pos[neighbors[j]]) - np.array(lcc_pos[node])

            # Calculate the angle between the two vectors
            angle = angle_between_vectors(vec1, vec2)

            # Store the angle as a node attribute (angles between edges meeting at this node)
            g.nodes[node][f'angle_{neighbors[i]}_{neighbors[j]}'] = angle

            # Store the angle between the two neighbors at the node in a dictionary
            lcc_angle_dict[(node, neighbors[i], neighbors[j])] = angle

            # Retrieve edge types for the two edges
            edge_1 = (node, neighbors[i]) if (node, neighbors[i]) in lcc_edge_types else (neighbors[i], node)
            edge_2 = (node, neighbors[j]) if (node, neighbors[j]) in lcc_edge_types else (neighbors[j], node)

            type_1 = lcc_edge_types[edge_1]
            type_2 = lcc_edge_types[edge_2]

            # Increment angle_id
            lcc_angle_id_counter += 1

            # Store the information in a list of tuples
            lcc_angle_data.append((lcc_angle_id_counter, edge_1, type_1, edge_2, type_2, angle))

# Create a DataFrame from the angle data
lcc_angle_df = pd.DataFrame(lcc_angle_data, columns=["angle_id", "edge_1", "type_1", "edge_2", "type_2", "angle_degree"])



################################
""" CALCULATE BOND STIFFNESS """
################################

# ---------------------------
# Calculate k_stretch from analytic Morse curvature (per-edge or scalar)
# ---------------------------
def k_from_morse_analytic(D, alpha):
    """
    analytic small-oscillation stiffness for Morse: k = 2 * D * alpha^2
    Inputs:
      D: well depth in energy units (same units as your sim energy, e.g., kB*T or J)
      alpha: Morse potential's "inverse-length parameter" (units 1/length in same length units as pos)
             ...we rescaled 3/alpha to match depletion potential, but units are still 1/length
             If alpha is scalar -> returns scalar; if array -> returns per-edge array
    Returns:
      k (same energy/length^2 units as D/length^2)
    """
    D = np.asarray(D, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    return 2.0 * D * (alpha**2)

# ---------------------------
# Calculate k_bend from r_ik
# ---------------------------
def compute_rik(angle_df, edge_df):
    """
    angle_df: DataFrame with columns:
        (node=j, neighbors i,k), angle_degree
    edge_df: DataFrame with columns:
        source, target, length (center–center length)

    Returns:
        angle_df with new columns: r_ik and kappa_perp
    """

    # Create a lookup for edge lengths
    edge_lookup = {}
    for _, row in edge_df.iterrows():
        a, b = row['source'], row['target']
        edge_lookup[(a, b)] = row['edge_length']
        edge_lookup[(b, a)] = row['edge_length']


    r_ik_list = []

    for _, row in angle_df.iterrows():
        edge_1, edge_2 = row['edge_1'], row['edge_2'],
        theta = np.deg2rad(row['angle_degree'])

        # lengths a = r_ij, b = r_jk
        a = edge_lookup[edge_1]
        b = edge_lookup[edge_2]

        # law of cosines
        rik = np.sqrt(a*a + b*b - 2*a*b*np.cos(theta))

        r_ik_list.append(rik)

    angle_df['r_ik'] = r_ik_list
    return angle_df


# calculate stretch stiffness and bend stiffness 
print(" - mean stretch and bend stiffness, and estiamted G':")

# If your network is not fully connected and you want to compare the results
# from the full network and the LCC, use "edge_df" and "angle_df" instead,
# to calculate stiffness from the full network

print("   LCC stiffness")
avg_rij = np.mean(edge_df_lcc['edge_length'].to_numpy())
avg_rij_sq = np.mean((edge_df_lcc['edge_length'].to_numpy())**2)
k_equip_r0 = round(kT/(avg_rij_sq),3)
var_rij = avg_rij_sq - avg_rij**2
k_equip_var_rij = round(kT/(var_rij),3)
print(f"   - k_stretch")
print(f"       - equipartition:")
print(f"         <r_ij> = {avg_rij}, <r_ij^2> {avg_rij_sq}, var_rij = {var_rij}")
print(f"         kT/(<r_ij^2>-<r_ij>^2) = {k_equip_var_rij}")

# compare to the analytical calculation of stiffness from the
# second derivative of the interaction potential
D_sim = U0 * kT #D0 * kT 
alpha_sim = kappa_m # DISTANCE
R_C_sim = 1 
r0_sim = 2*(R_C_sim+(rmin/550))
# analytic curvature:
k_morse_analytic = k_from_morse_analytic(D_sim, alpha_sim)
print("        - analytic morse curvature:", k_morse_analytic)

# calculate the bending stiffness from geometry
K_bend = m_K # units of [sim_energy] not units of [kT]
compute_rik(lcc_angle_df, edge_df_lcc)
mean_rik = np.mean(lcc_angle_df['r_ik'].to_numpy())
mean_rik_sq = np.mean(lcc_angle_df['r_ik'].to_numpy()**2)
var_rik = mean_rik_sq - mean_rik**2
k_bend_rik0 = (m_K)/(mean_rik_sq)
k_bend_var = (m_K)/(var_rik)
print(f"   - k_bend = {K_bend/kT}kT/<r_ik^2>")
print(f"          <r_ik> = {mean_rik}, <r_ik^2> = {mean_rik_sq}, var_rik = {var_rik}")
print(f"          K/(<r_ik^2>-<r_ik>^2)  = {k_bend_var}")

lcc_df['k_equip_var_rij'] = k_equip_var_rij
lcc_df['k_morse_analytic'] = k_morse_analytic 
lcc_df['k_bend_var'] = k_bend_var 

print(' - Save data:')

lcc_df.to_csv(data_outpath+'/lcc.csv',index = False)
print('   - Basic LCC network information calculated, saved to "lcc.csv"')



print(' - Angle analysis:')

lcc_angle_df.to_csv(data_outpath+'/angle_dist_lcc.csv',index = False)
angle_df.to_csv(data_outpath+'/angle_dist_all.csv',index = False)

tags = ['Full network:', 'LCC:']
angle_dfs = [angle_df, lcc_angle_df]
messages = [
    '      Full angle distribution saved to "angle_dist.csv"',
    '      LCC angle distribution saved to "lcc_angle_dist.csv"',
    ]

for a in range(len(angle_dfs)):

    tag = tags[a]
    df = angle_dfs[a]
    message = messages[a]
    print('    '+tag)
    ## OPTIONAL: plot angle distribution
    bond_types = list(pd.unique(edge_df['edge_type']))
    plot_message = ''

    nbins = 100
    # Plot the distribution of angles
    plt.hist(df['angle_degree'], bins=nbins, color='blue', alpha=0.7, density=True)

    # Convert the y-axis values to percentages
    plt.gca().yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1))

    plt.title('Distribution of Angles Between Edges')
    plt.xlabel('Angle (degrees)')
    plt.ylabel('Probability (%)')
    plt.ylim(top=0.03)
    plt.axvline(x=60, color='black', ls='--', lw=0.75)
    #plt.axvline(x=88.85,color='black',ls='--',lw=0.75)
    if 'LCC' in tag:
        plt.savefig(data_outpath+'/lcc_angledist.png', format='png', dpi=300, bbox_inches='tight')
        plot_message = ' and "lcc_angledist.png"'
    else:
        plt.savefig(data_outpath+'/angledist.png', format='png', dpi=300, bbox_inches='tight')
        plot_message = ' and "angledist.png"'
