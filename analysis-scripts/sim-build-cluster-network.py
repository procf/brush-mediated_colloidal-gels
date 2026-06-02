## Perform GMM clustering on a particle network 
## NOTE: expects data created by sim-build-particle-network.py  
## NOTE: this code assumes 1 colloid type (typeid=0)
##
"""
## Cluster the network using GMM
##
## How to use:
##     python sim-build-cluster-network.py brush_density
##
## INPUT:
## use lcc network (this is the elastic part) 
##    - node_df_lcc 
##        # keys : 'Tag#', 'x', 'y', 'z', 'Radius', 'TypeID'  
##        saved to data_outpath+'/node_df_lcc.csv' 
##    - edge_df_lcc 
##        # keys : 'source', 'target', 'edge_type'  
##        saved to data_outpath+'/edgelist_lcc.csv'
##    - g 
##        # keys : 'pos' [x,y,z], 'radius', 'type_id', 'edge_type' 
##        saved to data_outpath+'/lcc-particle-network_g.pkl'
##    - lcc_df 
##        # keys : 'cut_off', 'n_components', 'lcc_size', 'ncolloids', 'avg_degree', 
##        saved to data_outpath+'/lcc.csv'
##    - angle_df_lcc
##        # keys : 'angle_id', 'edge_1', 'type_1', 'edge_2', 'type_2', 'angle_degree' 
##        saved to data_outpath+'/angle_dist_lcc.csv'
##
## OUTPUT:
##    - embeddings10d 
##        saved to data_outpath+'/lcc_node_embeddings_10D.csv' (requires mapper to reload as a graph)
##    - bic_df 
##        # keys : 'k' 'BIC'  
##        saved to data_outpath+'/bic_scores.csv'
##    - clusters_df 
##        # keys : 'node', 'cluster' 
##        saved to 'clusters.csv'
##    - weighted_cluster_edges_df 
##        # keys : 'source_cluster' 'target_cluster' 'weight' 
##        saved to data_outpath+'/weighted_cluster_edges.csv'
##    - clustered_df, lcc_clustered_df
##        # keys : 'Tag#', 'x', 'y', 'z', 'Radius', 'TypeID', 'Cluster'
##        saved to data_outpath+'/clustered_df.csv'
##        saved to data_outpath+'/lcc_clustered_df.csv'
##    - cluster_diameters_df
##        # keys : 'Cluster' 'Diameter' 'Physical Diameter'
##        saved to data_outpath+'/cluster_diameters.csv'
##    - G_clusters
##        # clustered LCC network : 'cluster_id', 'size', 'diameter', 'physical_diameter'
##        saved to data_outpath+'/lcc-cluster-network_G_clusters.pkl'
##    - cluster_cauchyborn_df
##        # keys : 'average_xi', 'system_volume', 'phi_C', 'z_c'
##        saved to data_outpath+'/cluster_cauchy_born.csv'
"""
## (Rob Campbell)


########################
""" MODULE LIBRARY """
########################
# load data and build a network
import numpy as np
import pandas as pd
import networkx as nx
import os
import sys
import pickle

# for mapping bond information into network space
from node2vec import Node2Vec
import umap

# for clustering
from sklearn.mixture import GaussianMixture
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import ConvexHull
from scipy.spatial import cKDTree

# for faster BIC minimization
import multiprocessing
from multiprocessing import Pool

# for cluster analysis
from collections import defaultdict

# for plotting
import matplotlib.pyplot as plt



##########################
""" INPUT PARAMETERS """
##########################

brush_density = str(sys.argv[1]) # 'H' high 'L' low

if brush_density == 'H':
  gel_file = '../sim-scripts/potential-morsebrush/K_1xUC/hb_20nm-H/6mgmL-10mM-K10/Gelation-6mgmL-10mM-K10_H-BD.gsd'
  tag ='6mgmL-10mM-H_sim'
  hb = 20
  kappa_m = 75
  m_K = 10
  U0 = 9.7

if brush_density == 'L':
  gel_file = '../sim-scripts/potential-morsebrush/K_1xUC/hb_20nm-L/6mgmL-10mM-K10/Gelation-6mgmL-10mM-K10_L-BD.gsd'
  tag ='6mgmL-10mM-L_sim'
  hb = 20
  kappa_m = 75
  m_K = 20
  U0 = 20.4

# real parameters
kT = 0.1
d_g = 35 #[nm] = 2*r_g = 2*17.5
R_C_real = 550 #[nm]

# sim parameters
R_C = 1
PBC = True
cut_off = round((d_g+(2*hb))/R_C_real,3)
colloid_typeid = 0

# filepath to folder where data files will be created
data_outpath = 'data_'+tag

# create "data" subfolder if it doesn't exit
if os.path.exists(data_outpath) == False:
  os.mkdir(data_outpath)


#################
""" LOAD DATA """
#################

edge_df_lcc = pd.read_csv(data_outpath+'/edgelist_lcc.csv')
node_df_all = pd.read_csv(data_outpath+'/node_df_all.csv')
node_df_lcc = pd.read_csv(data_outpath+'/node_df_lcc.csv')
pos_lcc = node_df_lcc[['Tag#','x','y','z']]
lcc_df = pd.read_csv(data_outpath+'/lcc.csv')
ncolloids_lcc = lcc_df['ncolloids'].values[0]
L_X = lcc_df['L_X'].values[0]
L_Y = lcc_df['L_Y'].values[0]
L_Z = lcc_df['L_Z'].values[0]
Lbox = np.array([L_X, L_Y, L_Z])
with open(data_outpath + '/lcc-particle-network_g.pkl', 'rb') as fr:
    g = pickle.load(fr)

#####################
""" VECTORIZATION """
#####################
print(' - Starting GMM vectorization...')

## STEP 1: vectorize edgelist (convert to network space)
##   - takes ~1min for 1230 colloids
##   - saves data to CSV

# ensure node labels are integers
g = nx.relabel_nodes(g, lambda x: int(x))

embeddings_file_10D = data_outpath+'/lcc_node_embeddings_10D.csv'
if os.path.exists(embeddings_file_10D) == False:

  # (a)  convert from 3D real space to 128D network space with node2vec
  n_dims=128 # 128 dimensions to match Nabizadeh 2023
  node2vec = Node2Vec(g, dimensions=n_dims, walk_length=80, num_walks=20, workers=1)
  model = node2vec.fit(window=10, min_count=1)
  node_embeddings = model.wv
  print('   - '+str(n_dims)+'D node embeddings created')

  # [optional] save the full embeddings
  #embeddings_file_128D = data_outpath+'/node_embeddings.txt'
  #node_embeddings.save_word2vec_format(embeddings_file_128D)

  # (b)  dimensional reduction with UMAP (convert 128D to 10D)
  embeddings_128d = np.array([
      node_embeddings[str(int(node))] 
      for node in g.nodes() 
      if str(int(node)) in node_embeddings
  ])

  n_components = 10
  n_neighbors = 12
  min_dist = 0
  reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=n_components)
  embeddings_10d = reducer.fit_transform(embeddings_128d)
  print('   - embeddings converted from '+str(n_dims)+'D to 10D')

  # save the 10D embeddings
  embeddings_file_10D = data_outpath+'/lcc_node_embeddings_10D.csv'
  f = open(embeddings_file_10D,'w')
  j=0
  for item in g.nodes():
      f.write("%d,"%item)
      for i in range(n_components-1):
          f.write("%f,"%embeddings_10d[j,i])
      f.write("%f\n"%embeddings_10d[j,n_components-1])
      j+=1
  f.close()
  print('   - EMBEDDINGS SAVED TO:',embeddings_file_10D)

else:
    print('   - WARNING: 10D embeddings already exist, not re-embedding in 128D and re-performing dimensional reduction to 10D for this data.')


####################
""" MINIMIZE BIC """
####################
## STEP 2: minimize BIC to find optimal number of clusters
## takes 10-15min for 1000 k_values and 1230 colloids
print('   - Minimizing BIC (may take 20+min)...')

# A truly rigorous search tests all possibilities from 1 cluster to separated particles: k = range(1,ncolloids_lcc)
#   but most systems of our size have k_optimal < 1000 ; plot BIC curve to confirm minimum was found
min_nodes_per_cluster = 3
k_max = min(ncolloids_lcc//min_nodes_per_cluster, min(1000, ncolloids_lcc))  # capped for efficiency 
k_values = range(1, k_max)


# load embeddings
embeddings_file_10D = data_outpath+'/lcc_node_embeddings_10D.csv'
mapper = np.genfromtxt(embeddings_file_10D, delimiter=',')
embeddings_10d = mapper[:,1:]

bic_file = data_outpath+'/bic_scores.csv'
if os.path.exists(bic_file) == False:
    # Fit GMM for each k and calculate BIC
    bic_scores = []
    for k in k_values:
        gmm = GaussianMixture(n_components=k, random_state=42)
        gmm.fit(embeddings_10d)
        bic_scores.append(gmm.bic(embeddings_10d))

    # [optional] save BIC scores to CSV
    bic_df = pd.DataFrame({
        'k': list(k_values),
        'BIC': bic_scores
    })
    bic_df.to_csv(bic_file, index=False)
    print('   - BIC VALUES SAVED TO:',bic_file)

else:
    print('   - WARNING: BIC scores already exist, not re-minimizing BIC for this data.')


# Load BIC data to find the optimal k 
# (the number of clusters that give with the minimum BIC score)
bic_file = data_outpath+'/bic_scores.csv'
bic_results_df = pd.read_csv(bic_file)
bic_scores = list(bic_results_df['BIC'])
optimal_k = k_values[np.argmin(bic_scores)]
print(f"   - Optimal number of clusters (k) for 10D embeddings: {optimal_k}")

## [OPTIONAL] Plot the BIC scores
plt.figure(figsize=(8, 6))
plt.plot(bic_results_df['k'], bic_results_df['BIC'], marker='o', label='embeddings (10D)')
plt.xlabel('Number of clusters (k)')
plt.ylabel('BIC Score')
plt.title('BIC Score vs. Number of Clusters')
plt.axvline(x=optimal_k, c='royalblue')
plt.text(x=optimal_k+10, y=10000, s='k='+str(optimal_k), c='royalblue')
plt.legend()
plt.savefig(data_outpath+'/BIC_minimization.pdf', format='pdf', dpi=300, bbox_inches='tight')


#####################################
""" GMM CLUSTERING WITH OPTIMAL_K """
#####################################
## STEP 3: generate the cluster network from optimal_k

# load embeddings
embeddings_file_10D = data_outpath+'/lcc_node_embeddings_10D.csv'
mapper = np.genfromtxt(embeddings_file_10D, delimiter=',')
embeddings_10d = mapper[:,1:]

# Fit GMM with the optimal number of clusters
gmm_final = GaussianMixture(n_components=optimal_k, random_state=0, 
                            covariance_type='full', n_init=10).fit(embeddings_10d)

# Predict the cluster for each node
# `cluster_labels` is an array of cluster labels
cluster_labels = gmm_final.predict(embeddings_10d)

# [optional] reload network information if you do not have it from above
#g = pickle.load(data_outpath+'/particle-network_G.pkl')
#g = nx.relabel_nodes(g, lambda x: int(x))

# match nodes to their assigned clusters and collect cluster information
c = 0                              # index for looping through cluster_labels
cluster_dict = defaultdict(dict)   # cluster dictionary #1 (index:'node','cluster')
cluster_c = 0                      # cluster index
Clusters = {}                      # cluster dictionary #2 (index:nodes)
Nodes= []                          # list to collect nodes
Cluster_index = []                 # list to collect cluster indices

for node in g.nodes():
    # loop through the GMM results with this index
    idx = cluster_labels[c]
    # collect node and cluster information
    cluster_dict[cluster_c]['node'] = node
    cluster_dict[cluster_c]['cluster'] = idx
    Nodes.append(node)
    Cluster_index.append(idx+1) # start cluster numbering from 1 not 0
    # group nodes into clusters
    if idx in Clusters:
        Clusters[idx].add(node)
    else:
        Clusters[idx] = set()
        Clusters[idx].add(node)
        pass
    c += 1
    cluster_c += 1
    
Final_Cluster_Dict = {'Cluster':Cluster_index,
            'Node':Nodes}
clusters_df = pd.DataFrame.from_dict(Final_Cluster_Dict)

# save clustered data
cluster_file = data_outpath+'/clusters.csv'
clusters_df.to_csv(cluster_file,index = False)

print('   - clustering complete |',str(optimal_k),'clusters (#1-'+str(optimal_k)+')')
#clusters_df.head()

# group nodes into clusters
clusters = clusters_df['Cluster'].to_numpy()
nodes = clusters_df['Node'].to_numpy()

# use ALL colloids, unclustered particles (colloids not in the LCC) will be assigned cluster 0
sorted_clusters = np.zeros(ncolloids_lcc)
sorted_clusters[nodes] = clusters

# save full clustered network data
clustered_df = node_df_all.copy()
clustered_df['Cluster'] = sorted_clusters
clustered_df.to_csv(data_outpath+'/clustered_df.csv',index = False)

# remove unclustered particles (cluster = 0)
sorted_clusters = sorted_clusters[sorted_clusters != 0]

lcc_clustered_df = node_df_lcc.copy()
lcc_clustered_df['Cluster'] = sorted_clusters
lcc_clustered_df.to_csv(data_outpath+'/lcc_clustered_df.csv',index = False)

#######
## STEP 4: Verify clustering AND labeling
# Check that clusters are actually connected on the graph
#   - disconnected clusters suggest either:
#      (a) labeling was incorrect
#      (b) clusters need to be manually cleaned (outliers reassigned correctly)
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


# choose how to compute cluster volumes later (options provided)
def cluster_max_distance(unwrapped_positions):
    """Return true max pairwise Euclidean distance of the points."""
    if len(unwrapped_positions) < 2:
        return 0.0
    # pdist returns pairwise distances; max is what we want
    return float(np.max(pdist(unwrapped_positions)))

def cluster_convex_hull_volume(unwrapped_positions):
    """Return convex hull volume of the cluster (requires at least 4 non-coplanar points)."""
    if len(unwrapped_positions) < 4:
        return 0.0
    try:
        hull = ConvexHull(unwrapped_positions)
        return float(hull.volume)
    except Exception:
        return 0.0

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


def find_nearest_cluster(component_nodes, df, current_cluster_id, cluster_trees, position_columns=('x','y','z')):
    """
    Fast KD-tree version: find nearest cluster to a disconnected component.
    Returns ID of cluster with smallest minimum distance.
    """
    # positions of component particles
    comp_xyz = df[df['Tag#'].isin(component_nodes)][list(position_columns)].values

    min_dist = np.inf
    nearest_cluster = None

    # loop through other clusters
    for cluster in df['Cluster'].unique():
        cluster = int(cluster)
        if cluster == current_cluster_id:
            continue

        # get the xyz positions of that cluster
        other_xyz = df[df['Cluster'] == cluster][list(position_columns)].values

        # get KD-tree and query nearest neighbor for every node in component
        tree = cluster_trees[cluster]
        dists, _ = tree.query(comp_xyz)

        # find the minimum over all component nodes
        d = dists.min()

        if d < min_dist:
            min_dist = d
            nearest_cluster = cluster

    return nearest_cluster

# Check clustering and collect cluster data
connected_clusters = 0
unconnected_clusters = 0
unconnected_cluster_list = []
notes_list = []

# Ensure cluster IDs are ints (avoid float .0 keys)
lcc_clustered_df['Cluster'] = lcc_clustered_df['Cluster'].astype(int)

cluster_nodes_unique = np.sort(lcc_clustered_df['Cluster'].unique())

# Build initial KD-trees (updated once per cluster)
def build_cluster_trees(df):
    return {
        cid: cKDTree(df[df['Cluster'] == cid][['x','y','z']].values)
        for cid in df['Cluster'].unique()
    }

cluster_trees = build_cluster_trees(lcc_clustered_df)

# get ready to add more clusters if needed during reclustering
next_cluster_id = lcc_clustered_df['Cluster'].max() + 1


# check for disconnected clusters and update cluster assignment as needed:
for cluster_id in cluster_nodes_unique:
    mask = lcc_clustered_df['Cluster'] == cluster_id
    cluster_nodes = lcc_clustered_df.loc[mask, 'Tag#'].values

    # subgraph for this cluster
    cluster_subgraph = g.subgraph(cluster_nodes)

    # ensure node labels are integers
    cluster_subgraph = nx.relabel_nodes(cluster_subgraph, lambda x: int(x))

    # skip fully connected clusters
    if nx.is_connected(cluster_subgraph):
        connected_clusters += 1
        continue

    # split or reassign disconnected components
    unconnected_clusters += 1
    notes = []

    components = list(nx.connected_components(cluster_subgraph))
    sizes = np.array([len(c) for c in components])
    total_size = sizes.sum()

    # find the LCC
    lcc_nodes = max(components, key=len)

    # get information about the disconnect
    lcc_size = len(lcc_nodes)
    n_disconnected = total_size - lcc_size
    percent_disconnected = round(n_disconnected / total_size * 100, 2)

    if percent_disconnected <= 10:
        correction = "reassign"

        # loop through each cluster cc
        for comp_nodes in components:
            comp_nodes = set(comp_nodes)

            # keep the LCC of this cluster unchanged
            if comp_nodes == lcc_nodes:
                continue

            # reassign small components to the nearest other cluster
            nearest = find_nearest_cluster(
                component_nodes=comp_nodes,
                df=lcc_clustered_df,
                current_cluster_id=cluster_id,
                cluster_trees=cluster_trees
            )
            notes.append(
                f"Moved {len(comp_nodes)}-particle fragment from Cluster {cluster_id} → {nearest}"
            )

            lcc_clustered_df.loc[
                lcc_clustered_df['Tag#'].isin(comp_nodes),
                'Cluster'
            ] = nearest

        # after all reassignments, update KD-trees 
        cluster_trees = build_cluster_trees(lcc_clustered_df)

    else:
        correction = "split"
        added_clusters = 0

        for comp_nodes, comp_size in zip(components, sizes):

            comp_nodes = set(comp_nodes)
            size_fraction = comp_size / total_size

            # keep the LCC as this cluster
            if comp_nodes == lcc_nodes:
                continue

            # save large components as new clusters with a new cluster ID
            if size_fraction > 0.10:

                lcc_clustered_df.loc[
                    lcc_clustered_df['Tag#'].isin(comp_nodes),
                    'Cluster'
                ] = next_cluster_id

                notes.append(
                    f"Created new cluster {next_cluster_id} from {len(comp_nodes)}-particle fragment of Cluster {cluster_id}"
                )

                next_cluster_id += 1
                added_clusters += 1

            # reassign small components to neighboring clusters
            else:
                nearest = find_nearest_cluster(
                    component_nodes=comp_nodes,
                    df=lcc_clustered_df,
                    current_cluster_id=cluster_id,
                    cluster_trees=cluster_trees
                )
                notes.append(
                    f"Reassigned tiny {len(comp_nodes)}-particle fragment of Cluster {cluster_id} → {nearest}"
                )

                lcc_clustered_df.loc[
                    lcc_clustered_df['Tag#'].isin(comp_nodes),
                    'Cluster'
                ] = nearest

        notes.append(f"Added {added_clusters} new clusters")
        # Update KD-trees after reassignments
        cluster_trees = build_cluster_trees(lcc_clustered_df)

    unconnected_cluster_list.append({
        'Cluster': cluster_id,
        'Size': len(cluster_subgraph),
        'LCC_Size': lcc_size, 
        'N_Disconnected': n_disconnected,
        '%_Disconnected': percent_disconnected,
        'correction': correction,
        'notes':notes
        })

if unconnected_clusters == 0:
    print(f"     {connected_clusters} Clusters (all fully-connected)")
else:
    percent_connection_error = round(unconnected_clusters/(unconnected_clusters+connected_clusters),2)*100
    print('     WARNING: disconnected clusters suggests an error occured...') 
    print(f"     {unconnected_clusters} Disconnected Clusters ({percent_connection_error}% clustering error)")

    # Convert the disconnected cluster data to a df for easy viewing
    unconnected_clusters_df = pd.DataFrame(unconnected_cluster_list)
    disconnected_data = unconnected_clusters_df[['Cluster','Size','LCC_Size','N_Disconnected','%_Disconnected','correction']]
    notes_array = unconnected_clusters_df['notes'].to_numpy()
    clusters = pd.unique(unconnected_clusters_df['Cluster'])
    for c in range(len(clusters)):
        cluster = clusters[c]
        print(disconnected_data.loc[disconnected_data['Cluster'] == cluster])
        notes_cluster = notes_array[c] 
        for note in notes_cluster:
            print(f"     -- {note}")

    disconnected_data.to_csv(data_outpath+'/unconnected_clusters.csv',index = False)

    # find new cluster info
    cluster_nodes_unique = np.unique(lcc_clustered_df['Cluster'])
    print(f"     ---")
    print(f"     {len(cluster_nodes_unique)} Clusters after processing (all fully-connected)")


physical_cluster_diameters = {}
cluster_diameters = {}
cluster_hull_volumes = {}
cluster_particle_counts = {}
cluster_inner_degree = {}
cluster_center_x = {}
cluster_center_y = {}
cluster_center_z = {}
cluster_angle_dfs = []

for cluster_id in cluster_nodes_unique:
    mask = lcc_clustered_df['Cluster'] == cluster_id
    cluster_nodes = lcc_clustered_df.loc[mask, 'Tag#'].values
    cluster_pos_df = lcc_clustered_df.loc[mask, ['Tag#','x','y','z']]

    # subgraph for this cluster
    cluster_subgraph = g.subgraph(cluster_nodes)

    # ensure node labels are integers
    cluster_subgraph = nx.relabel_nodes(cluster_subgraph, lambda x: int(x))

    topological_diameter = nx.diameter(cluster_subgraph)
    # use whole cluster for geometry if you prefer
    geometry_nodes = list(cluster_nodes)
    geometry_edges = list(cluster_subgraph.edges())

    # unwrap chosen geometry nodes
    if len(geometry_nodes) > 0:
        if PBC == True:
            unwrap_pos = unwrap_cluster_positions(cluster_pos_df, Lbox, geometry_edges, list(geometry_nodes))
            # compute true max distance (physical diameter = max pairwise distance + 2*particle_radius)
        if PBC == False:
            pos_df_indexed = cluster_pos_df.set_index("Tag#")
            unwrap_pos = pos_df_indexed.loc[list(geometry_nodes)][["x","y","z"]].values.astype(float)

        maxpair = cluster_max_distance(unwrap_pos)
        # if you assumed R_particle = 1, add 2*R (or use actual particle radii if available)
        # replace PARTICLE_RADIUS with your real radius if needed
        PARTICLE_RADIUS = 1.0
        physical_diameter = maxpair + 2.0 * PARTICLE_RADIUS

        # convex hull volume (optional)
        hull_vol = cluster_convex_hull_volume(unwrap_pos)

        # particle count
        n_particles = len(geometry_nodes)

        # cluster center (in unwrapped coordinates)
        cluster_center = unwrap_pos.mean(axis=0)

        if PBC == True:
          # wrap back into simulation box for placement of coarse-grained particle
          cluster_center_wrapped = cluster_center % Lbox
          cluster_center = cluster_center_wrapped.copy()

    else:
        physical_diameter = 0.0
        hull_vol = 0.0
        n_particles = 0
        cluster_center = np.array([0,0,0])

    # get interior average degree of each cluster
    nnodes_cluster = cluster_subgraph.number_of_nodes()
    nedges_cluster = nx.number_of_edges(cluster_subgraph)
    inner_avg_degree = 2 * nedges_cluster / nnodes_cluster

    # get interior angle information for each cluster
    # Prepare a list to store angle data
    angle_data = []

    # Retrieve node positions and edge types
    pos = nx.get_node_attributes(cluster_subgraph, 'pos')
    edge_types = nx.get_edge_attributes(cluster_subgraph, 'edge_type')

    # Angle ID counter
    angle_id_counter = 0

    # Calculate the angle between edges at each node
    angle_dict = {}
    for node in cluster_subgraph.nodes():
        neighbors = list(cluster_subgraph.neighbors(node))
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
                cluster_subgraph.nodes[node][f'angle_{neighbors[i]}_{neighbors[j]}'] = angle

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
    cluster_angle_dfs.append(angle_df)


    cluster_diameters[int(cluster_id)] = int(topological_diameter)
    physical_cluster_diameters[int(cluster_id)] = float(physical_diameter)
    cluster_hull_volumes[int(cluster_id)] = float(hull_vol)
    cluster_particle_counts[int(cluster_id)] = int(n_particles)
    cluster_inner_degree[int(cluster_id)] = float(inner_avg_degree)
    cluster_center_x[int(cluster_id)] = cluster_center[0]
    cluster_center_y[int(cluster_id)] = cluster_center[1]
    cluster_center_z[int(cluster_id)] = cluster_center[2]


# DataFrame
cluster_diameters_df = pd.DataFrame.from_dict(physical_cluster_diameters, orient='index', columns=['Physical Diameter'])
cluster_diameters_df['Cluster'] = cluster_diameters_df.index
cluster_diameters_df['Topological Diameter'] = pd.Series(cluster_diameters)
cluster_diameters_df['ConvexHullVolume'] = pd.Series(cluster_hull_volumes)
cluster_diameters_df['N_particles'] = pd.Series(cluster_particle_counts)
cluster_diameters_df['Inner_Avg_Degree'] = pd.Series(cluster_inner_degree)
cluster_diameters_df['CenterX'] = pd.Series(cluster_center_x)
cluster_diameters_df['CenterY'] = pd.Series(cluster_center_y)
cluster_diameters_df['CenterZ'] = pd.Series(cluster_center_z)
cluster_diameters_df.reset_index(drop=True, inplace=True)
cluster_diameters_df.to_csv(data_outpath + '/cluster_diameters.csv', index=False)




#############################
""" BUILD CLUSTER NETWORK """
#############################
## STEP 5: Create a new graph G_clusters with weighted edges

# set df to be indexed by Tag# for easier searching
lcc_clustered_df = lcc_clustered_df.set_index("Tag#")

# create a new graph of the clusters!
G_clusters = nx.Graph()

# Add nodes for each cluster scaled by cluster size (number of particles)
for cluster_id in np.unique(lcc_clustered_df['Cluster']):
    cluster_size = np.sum(lcc_clustered_df['Cluster'] == cluster_id)
    c_df = cluster_diameters_df.loc[cluster_diameters_df['Cluster'] == cluster_id]
    cluster_diameter = c_df['Topological Diameter'].values[0]
    physical_cluster_diameter = c_df['Physical Diameter'].values[0]
    pos = (c_df['CenterX'].values[0], c_df['CenterY'].values[0], c_df['CenterZ'].values[0])
    G_clusters.add_node(cluster_id, size=cluster_size, pos=pos, diameter=cluster_diameter, physical_diameter=physical_cluster_diameter)

# Add edges between clusters (this is the number of particle-particle contacts that connect each cluster)
for index, row in edge_df_lcc.iterrows():
    source_cluster = lcc_clustered_df.loc[row['source'], 'Cluster']
    target_cluster = lcc_clustered_df.loc[row['target'], 'Cluster']
    
    if source_cluster != target_cluster:
        if G_clusters.has_edge(source_cluster, target_cluster):
            G_clusters[source_cluster][target_cluster]['weight'] += 1
        else:
            G_clusters.add_edge(source_cluster, target_cluster, weight=1)

# Create a list to store the edges with their weights
edges_with_weights = []

# Iterate over the edges in G_clusters to extract source, target, and weight
for (source_cluster, target_cluster, weight) in G_clusters.edges(data='weight'):
    edges_with_weights.append({
        'source_cluster': int(source_cluster),
        'target_cluster': target_cluster,
        'weight': weight
    })

# Convert the list to a DataFrame
weighted_cluster_edges_df = pd.DataFrame(edges_with_weights)

if PBC == True:
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

    add_edge_lengths(G_clusters, box=Lbox) # with PBC
    # Extract all edge lengths from the Graph
    rows = []
    for u, v, data in G_clusters.edges(data=True):
        rows.append({
            "source_cluster": u,
            "target_cluster": v,
            "edge_length": data.get("length", np.nan)
        })

if PBC == False:
    def add_edge_lengths(G, delta_max, beta_U0):
        """
        Assign bond lengths from Boltzmann-weighted depletion shell.

        Parameters
        ----------
        G : networkx.Graph
        r_contact : float
            Contact separation
        delta_max : float
            Attraction range (Δ)
        beta_U0 : float
            Dimensionless well depth (β U0)
        """

        # Boltzmann slope
        k = beta_U0 / delta_max

        # normalization constant for truncated exponential
        norm = 1.0 - np.exp(-k * delta_max)

        for u, v, data in G.edges(data=True):

            # inverse transform sampling
            u_rand = np.random.rand()
            delta = -np.log(1 - u_rand * norm) / k

            R_u = G.nodes[u].get('R', 1.0) # defaults to 1 if there is no 'R'
            R_v = G.nodes[v].get('R', 1.0)

            data['length'] = (R_u + R_v) + delta ## center-center distance

    add_edge_lengths(G_clusters, 2*d_g*1e-3, U0) # with PBC
    # Extract all edge lengths from the Graph
    rows = []
    for u, v, data in G_clusters.edges(data=True):
        rows.append({
            "source_cluster": u,
            "target_cluster": v,
            "edge_length": data.get("length", np.nan)
        })

edge_len_df = pd.DataFrame(rows)
# Merge the lengths back into edge_df
weighted_cluster_edges_df = weighted_cluster_edges_df.merge(edge_len_df,
                       on=["source_cluster", "target_cluster"],
                       how="left")

weighted_cluster_edges_df.to_csv(data_outpath+'/weighted_cluster_edges.csv',index = False)

# Display the DataFrame
#print(weighted_cluster_edges_df.head())

# SAVE GRAPH AS PKL FILE
particle_network_filename = "lcc-cluster-network_G_clusters.pkl"
with open(data_outpath+'/'+particle_network_filename, 'wb') as fr:
  pickle.dump(G_clusters, fr)
print(' - Graph saved to "lcc-cluster-network_G_clusters.pkl" (PKL file)')


# cluster-angles
# Prepare a list to store angle data
cluster_angle_data = []

cluster_pos = nx.get_node_attributes(G_clusters, 'pos')

# Angle ID counter
cluster_angle_id_counter = 0

# ensure node labels are integers
G_clusters = nx.relabel_nodes(G_clusters, lambda x: int(x))

# Calculate the angle between edges at each node
cluster_angle_dict = {}
for node in G_clusters.nodes():
    neighbors = list(G_clusters.neighbors(node))
    if len(neighbors) < 2:
        continue  # No angle if less than two edges meet at the node

    for i in range(len(neighbors)):
        for j in range(i + 1, len(neighbors)):
            # Get position vectors for two neighbors
            if PBC == True:
              vec1 = pbc_vec(np.array(cluster_pos[node]), np.array(cluster_pos[neighbors[i]]), Lbox) # PBC
              vec2 = pbc_vec(np.array(cluster_pos[node]), np.array(cluster_pos[neighbors[j]]), Lbox) # PBC
            else:
              vec1 = np.array(cluster_pos[neighbors[i]]) - np.array(cluster_pos[node])
              vec2 = np.array(cluster_pos[neighbors[j]]) - np.array(cluster_pos[node])

            # Calculate the angle between the two vectors
            angle = angle_between_vectors(vec1, vec2)

            # Store the angle as a node attribute (angles between edges meeting at this node)
            G_clusters.nodes[node][f'angle_{neighbors[i]}_{neighbors[j]}'] = angle

            # Store the angle between the two neighbors at the node in a dictionary
            angle_dict[(node, neighbors[i], neighbors[j])] = angle

            # Increment angle_id
            cluster_angle_id_counter += 1

            # Retrieve edge types for the two edges
            edge_1 = (node, neighbors[i]) 
            edge_2 = (node, neighbors[j]) 

            # Store the information in a list of tuples
            cluster_angle_data.append((cluster_angle_id_counter, edge_1, edge_2, angle))

# Create a DataFrame from the angle data
cluster_angle_df = pd.DataFrame(cluster_angle_data, columns=["angle_id", "edge_1", "edge_2", "angle_degree"])

## bin cluster sizes to determine the cluster length scale xi

# average length scale
diams = cluster_diameters_df['Physical Diameter']
average_xi = sum(diams) / len(diams)
bin_size = 1 # bin physical size by integer distances

# Define the number of bins or specific bin edges
max_diam = int(np.ceil(cluster_diameters_df['Physical Diameter'].max()))
bins = np.linspace(0,max_diam,max_diam+1)

counts, bin_edges = np.histogram(
    cluster_diameters_df['Physical Diameter'],
    bins=bins,
    density=False
)

# Normalize to probability
normalized_counts = counts / counts.sum()

bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# Plot the results
plt.figure(figsize=(10, 6))
plt.bar(bin_centers, normalized_counts, width=(bin_edges[1]-bin_edges[0]))
plt.axvline(x=average_xi,c='black')
plt.text(x=average_xi+0.2,y=0.14,s='$\\langle \\xi \\rangle=$'+str(round(average_xi,2))+' edges', c='black', size=15)

plt.xlabel('$D_c$ (cluster diameter, [$D/D_{colloid}$])',fontsize=15)
plt.ylabel('$P(D_c)$',fontsize=15)
plt.title('Cluster Size Distribution',fontsize=15)
plt.savefig(data_outpath+'/cluster_sizedist.pdf', format='pdf', dpi=300, bbox_inches='tight')

diams = np.array(list(physical_cluster_diameters.values()))
average_xi = sum(diams) / len(diams)
#print(' - xi:',average_xi) # in colloid radii

system_volume = Lbox[0] * Lbox[1] * Lbox[2]

# cluster volume fraction
# volume of particles
R_p = 1
particle_vol = 4.0/3.0 * np.pi * R_p**3
total_cluster_particle_volume = sum(cluster_particle_counts.values()) * particle_vol
system_volume = Lbox[0] * Lbox[1] * Lbox[2]
phi_C_particles = total_cluster_particle_volume / system_volume
#print(' - phi_C_particles:',phi_C_particles)

# volume from convex hull
phi_C_hull_sum = sum(cluster_hull_volumes.values()) / system_volume
#print(' - phi_C_hull_sum:',phi_C_hull_sum)

# as spheres
diams = np.array(list(physical_cluster_diameters.values()))
volumes = 4.0/3.0 * np.pi * (0.5 * diams)**3
phi_C_spheres = np.sum(volumes) / system_volume
#print(' - phi_C:',phi_C_spheres)
#phi_C = sum(4*np.pi*(0.5*np.array(diams))**3/3)/system_volume
#print(' - phi_C:',phi_C)

#print("CLUSTER VOLUME FRACTIONS")
#print((np.array(cluster_particle_counts.values()) * particle_vol)/volumes)


# cluster average coordination number
n_clusters = G_clusters.number_of_nodes()
if n_clusters > 0:
    total_edges = G_clusters.number_of_edges()   # if graph is simple, E
    avg_degree_unweighted = 2.0 * total_edges / n_clusters

    # if you have weights and want weighted average degree:
    weighted_degrees = dict(G_clusters.degree(weight='weight'))
    avg_degree_weighted = np.mean(list(weighted_degrees.values()))
else:
    avg_degree_unweighted = 0.0
    avg_degree_weighted = 0.0

#print(' - z_c (unweighted avg degree):', avg_degree_unweighted)
#print(' - z_c (weighted avg degree):', avg_degree_weighted)

# compile outputs
data={
      'average_xi'             :[average_xi],
      'system_volume'          :[system_volume],
      'phi_C_particles'        :[phi_C_particles],
      'phi_C_hull_sum'         :[phi_C_hull_sum],
      'phi_C_spheres'          :[phi_C_spheres],
      'z_c_unique'             :[avg_degree_unweighted],
      'z_c_weighted'           :[avg_degree_weighted],
      }

cluster_cauchyborn_df = pd.DataFrame(data)

# STEP 1: ESTIMATE CALLADINE RELATION
def compute_calladine_F_over_N(z_mean, d=3, r=None, c=None, extra_constraints_per_particle=0.0,
                              N_triv=0, N_particles=None):
    """
    Compute generalized (F - S)/N ~ (d + r) - (z_mean * c)/2 - extra_constraints_per_particle - N_triv/N.
    """
    if r is None:
        r = d * (d - 1) // 2
    if c is None:
        c = 1
    base = (d + r) - 0.5 * z_mean * c - extra_constraints_per_particle
    if (N_particles is not None) and (N_particles > 0):
        base -= float(N_triv) / float(N_particles)
    return base


clusters = cluster_diameters_df['Cluster']
inner_data_dfs = []
for c in range(len(clusters)):
    cluster = clusters[c]
    avg_degree = cluster_diameters_df.loc[cluster_diameters_df['Cluster'] == cluster]['Inner_Avg_Degree'].values

    angle_df = cluster_angle_dfs[c]
    n_angles = len(angle_df)

    m = n_angles / cluster_diameters_df.loc[cluster_diameters_df['Cluster'] == cluster]['N_particles'].values[0]

    # 3D real-space system
    N_triv_free = 6
    #periodic system or LARGE system
    N_triv_free = 0

    F_over_N_frictionless = compute_calladine_F_over_N(
          avg_degree, d=3, r=0, c=1, extra_constraints_per_particle=0,
          N_triv=N_triv_free, N_particles=ncolloids_lcc
    )
    F_over_N_frictional = compute_calladine_F_over_N(
          avg_degree, d=3, r=3, c=3, extra_constraints_per_particle=0,
          N_triv=N_triv_free, N_particles=ncolloids_lcc
    )

    # angle constraints
    F_over_N_bending = compute_calladine_F_over_N(
          #avg_degree, d=3, r=3, c=3, extra_constraints_per_particle=1.5,
          avg_degree, d=3, r=0, c=1, extra_constraints_per_particle=m,
          N_triv=N_triv_free, N_particles=ncolloids_lcc
    )

    df = pd.DataFrame([n_angles], columns=['n_angles'])
    df['m'] = m
    df['MC_frictionless'] = F_over_N_frictionless
    df['MC_frictional'] = F_over_N_frictional
    df['MC_bending'] = F_over_N_bending
    df['z_iso'] = 6-2*m
    inner_data_dfs.append(df)

    #print(f" - Cluster {cluster}: <z>={avg_degree[0]}, n_angles={n_angles}")
    #print(f"     - (F-S)/N = {F_over_N_frictionless} (frictionless)")
    #print(f"     - (F-S)/N = {F_over_N_frictional} (frictional)")
    #print(f"     - (F-S)/N = {F_over_N_bending} (bending)")

cluster_interior_df = pd.concat(inner_data_dfs, ignore_index=True)
#print(cluster_interior_df)
cluster_diameters_df = pd.concat([cluster_diameters_df, cluster_interior_df], axis=1)
cluster_diameters_df.to_csv(data_outpath + '/cluster_diameters.csv', index=False)

## Load k_stretch and k_bend 
k_stretch = lcc_df['k_equip_var_rij'].values[0]
k_bend = lcc_df['k_bend_var'].values[0]
print(f'k_stretch from particles = {k_stretch}')
print(f'k_bend from particles = {k_bend}')


###############################################################
# FILTERED FOR RIGIDITY — clusters with Inner_Avg_Degree >= 2.4
###############################################################

# 1. Identify clusters that satisfy the threshold
threshold = 2.4

# If you have the info in a column of cluster_diameters_df:
filtered_clusters = cluster_diameters_df.loc[
    cluster_diameters_df["Inner_Avg_Degree"] >= threshold, "Cluster"
].unique()


# Safety check
if len(filtered_clusters) == 0:
    print("No clusters satisfy Inner_Avg_Degree >= 2.4")
    filtered_stats = {
        "average_xi_filtered": np.nan,
        "phi_C_spheres_filtered": np.nan,
        "z_c_unique_filtered": np.nan,
        "z_c_weighted_filtered": np.nan
    }
else:

    # average_xi_filtered
    diams_filtered = np.array([
        physical_cluster_diameters[cid]
        for cid in filtered_clusters
        if cid in physical_cluster_diameters
    ])
    average_xi_filtered = diams_filtered.mean()

    # phi_C_spheres_filtered (sum volume of selected clusters)
    volumes_filtered = 4.0/3.0 * np.pi * (0.5 * diams_filtered)**3
    phi_C_spheres_filtered = np.sum(volumes_filtered) / system_volume

    # Build subgraph containing only filtered clusters
    G_clusters_filtered = G_clusters.subgraph(filtered_clusters).copy()

    # Compute z_c_unique_filtered and z_c_weighted_filtered
    n_filt = G_clusters_filtered.number_of_nodes()
    if n_filt > 0:
        E_filt = G_clusters_filtered.number_of_edges()

        # Unweighted mean degree
        z_c_unique_filtered = 2.0 * E_filt / n_filt

        # Weighted mean degree
        weighted_deg_filt = dict(G_clusters_filtered.degree(weight='weight'))
        z_c_weighted_filtered = np.mean(list(weighted_deg_filt.values()))
    else:
        z_c_unique_filtered = 0.0
        z_c_weighted_filtered = 0.0

    # results
    filtered_stats = {
        "average_xi": average_xi_filtered,
        "phi_C_spheres": phi_C_spheres_filtered,
        "z_c_unique": z_c_unique_filtered,
        "z_c_weighted": z_c_weighted_filtered
    }

# use the data that has been filtered for rigid clusters
iso_df = pd.DataFrame([filtered_stats])


# xi = cluster diameter
G_c = (4*k_stretch*iso_df['z_c_weighted'].values[0]*iso_df['phi_C_spheres'].values[0])/(5*np.pi*iso_df['average_xi'].values[0])
G_b = (124*k_bend*iso_df['z_c_weighted'].values[0]*iso_df['phi_C_spheres'].values[0])/(135*np.pi*iso_df['average_xi'].values[0])
print(f"   - Cauchy-Born (z_c == all particle contacts between clusters)")
print(f"     SIM UNITS")
print(f"     - G' = {round((G_c+G_b),2)}")
print(f"     - G'_c = {round(G_c,2)}")
print(f"     - G'_b = {round(G_b,2)}")

iso_df['G_prime'] = round(G_c+G_b,2)
iso_df['G_c'] = G_c
iso_df['G_b'] = G_b

print(f"z_c_unique = {iso_df['z_c_unique'].values[0]}")
print(f"z_c_weighted = {iso_df['z_c_weighted'].values[0]}")
print(f"phi_C = {iso_df['phi_C_spheres'].values[0]}")
print(f"xi = {iso_df['average_xi'].values[0]}")

iso_df.to_csv(data_outpath+'/cluster_cauchy_born.csv',index = False)
