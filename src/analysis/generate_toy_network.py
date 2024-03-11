import os
import argparse
import pickle as pk
import sys
import re

import networkx as nx
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

sys.path.append('src/utils')
from loaders import load_signaling_df

def gc_id_to_resid(gc_id: str) -> int:
    """Gets a residue id from the getcontacts frequency file and returns only the residue index

    Args:
        gc_id (str): residue identifier consisting in <chain>:<resname>:<resid>

    Returns:
        resid (int): Residue index of the residue
    """
    assert type(gc_id) is str
    # Generate a list from the getcontacts id separating the fields by ":"
    gc_id_list = gc_id.split(':')
    assert len(gc_id_list) == 3
    # Extract the residue id from the last element of the list identifier
    resid = gc_id_list[-1]
    assert resid.isnumeric()
    resid = int(resid)
    return resid

def get_graph_from_freqs(path):
    # Load interaction frequencies from wt
    freq_df = pd.read_csv(path, sep='\t', skiprows=2, header=None)
    freq_df.columns = ['r1', 'r2', 'freq']
    freq_df.r1 = freq_df.r1.apply(gc_id_to_resid)
    freq_df.r2 = freq_df.r2.apply(gc_id_to_resid)
    # Filter interactions from conscutive residues
    no_consecutives_mask = ~((freq_df.r1 - freq_df.r2).abs() == 1)
    freq_df = freq_df.loc[no_consecutives_mask]
    g = nx.graph.Graph()
    freq_df.loc[:, 'weight'] = 1/freq_df.freq.values
    edges = []

    for _, (r1, r2, _, weight) in freq_df.iterrows():
        edges.append((int(r1), int(r2), {'weight': weight}))
    g.add_edges_from(edges)
    
    return g

def process_resid(res_id):
    residx = int(res_id.split(':')[-1])
    return residx

def parse_pathway_data(path):
    
    df = pd.read_csv(path, index_col=0)
    df = df.reset_index()
    df = df.melt(id_vars='index', var_name='r2')
    df = df.rename(columns={'index': 'r1'})
    df = df[df.value != 0]
    df.r1 = df.r1.apply(process_resid)
    df.r2 = df.r2.apply(process_resid)
    return df

def preprocess_graph(graph_path):
    """Takes a path to a pickled graph and returns a graph with reformated node names"""

    # Assert that all nodes are strings of the format <chain>:<resname>:<resid>
    G = pk.load(open(graph_path, 'rb'))
    for node in G.nodes:
        assert type(node) is str, "Node names are not strings"
        # assert re.match(r'[A-Z]:[A-Z]{3}:\d+', node), f"Node names are not in the format chain:resname:resid ({node})"
    # Rename nodes to only contain the residue index
    mapping = {old_name: process_resid(old_name) for old_name in G.nodes}
    G = nx.relabel_nodes(G, mapping, copy=True)
    return G

def path_df_to_list(path_df, threshold):
    """Takes a dataframe with the ACN and returns a list of nodes in the path"""
    
    filgered_df = path_df[path_df.value > threshold]
    shortest_path_nodes = np.unique(filgered_df[['r1', 'r2']].values.flatten())
    shortest_path_nodes = shortest_path_nodes.tolist()
    return shortest_path_nodes

def get_gi_preferred_mutants():
    signaling_df = load_signaling_df()
    gi_preferred_mutants = signaling_df.index[signaling_df.profile == 1].astype(int).values
    return gi_preferred_mutants
    

def get_distances_to_shortest_pathway(g: nx.graph.Graph,
                                      shortest_pathway: list,
                                      node_set: list,
                                      weighted=False) -> pd.DataFrame:
    
    # Validate input
    # Check that all nodes in shortest_pathway and node_list are in g
    assert all([node in g for node in shortest_pathway])
    # Check that node_lists have lengths different than 0
    assert len(shortest_pathway) > 0
    assert len(node_set) > 0

    # Create list of shortest path lengths
    shortest_path_lengths = []
    for node in node_set:
        if node not in g:
            # this can happen if the node corresponds to a residue that has not been crystallized
            continue
        # Create list for storing lengths
        node_shortest_path_lengths = []
        if weighted:
            shortest_path_length = nx.shortest_path_length(g, node, weight='weight')
        else:
            shortest_path_length = nx.shortest_path_length(g, node)
        
        for path_node in shortest_pathway:
            # compute shortest path length
            if path_node not in shortest_path_length:
                continue
            node_shortest_path_lengths.append(shortest_path_length[path_node])
        # Select lowest path length
        shortest_path_lengths.append([min(node_shortest_path_lengths), node])
        
    shortest_path_lengths = pd.DataFrame(shortest_path_lengths, columns=['distance', 'residue'])
        
    return shortest_path_lengths
            
def gen_random_data(graph: nx.graph.Graph, n_nodes: int, shortest_path: list, n_replicates=100, weighted=False):
    
    all_shortest_distances = []
    
    for replica in range(n_replicates):
        
        # Select n_nodes random_nodes from the graph
        random_nodes = np.random.choice(graph.nodes, n_nodes)
        
        shortest_distance_list = get_distances_to_shortest_pathway(graph, shortest_path, random_nodes, weighted=weighted)['distance']
        
        all_shortest_distances += list(zip(shortest_distance_list, [replica]*len(shortest_distance_list)))
        
    all_shortest_distances_df = pd.DataFrame(all_shortest_distances, columns=['distance', 'replica'])
    
    return all_shortest_distances_df
        

    
def main(graph_path, shortest_path_path, threshold, replicate):
    
    threshold_str = str(threshold).split('.')[1][:3]
    # extract the acn results path 
    results_dir = '/'.join(shortest_path_path.split('/')[:-2])
    df_path = f'{results_dir}/allosteric_network_distance/dist_to_pathway_r{replicate}_t{threshold_str}.csv'
    random_df_path = f'{results_dir}/allosteric_network_distance/random_dist_to_pathway_r{replicate}_t{threshold_str}.csv'

    graph = preprocess_graph(graph_path)
    
    shortest_pathway_df = parse_pathway_data(shortest_path_path)
    shortest_pathway = path_df_to_list(shortest_pathway_df, float(threshold))
    

    gi_preferred_mutants = get_gi_preferred_mutants()
    
    weighted_df = get_distances_to_shortest_pathway(graph, shortest_pathway, gi_preferred_mutants, weighted=True)
    df = get_distances_to_shortest_pathway(graph, shortest_pathway, gi_preferred_mutants, weighted=False)

    # Generate Null hypothesis random data
    n_nodes = len(gi_preferred_mutants)
    random_data_df = gen_random_data(graph, n_nodes, shortest_pathway, n_replicates=1000, weighted=False)
    random_data_df['type'] = 'unweighted'
    weighted_random_data_df = gen_random_data(graph, n_nodes, shortest_pathway, n_replicates=1000, weighted=True)
    weighted_random_data_df['type'] = 'weighted'
    random_data_df = pd.concat([random_data_df, weighted_random_data_df])
    
    weighted_df['theshold'] = threshold
    weighted_df['type'] = 'weighted'

    df['theshold'] = threshold
    df['type'] = 'unweighted'
    
    df = pd.concat([weighted_df, df])

    os.makedirs(f'{results_dir}/allosteric_network_distance', exist_ok=True)
    
    df.to_csv(df_path)
    
    random_data_df.to_csv(random_df_path)
    
    # # Plot data
    # plt.figure(dpi=300)
    # sns.histplot(random_data, label='Random distribution')
    # # sns.histplot(shorteset_path_lengths, label='Distances')
    # plt.axvline(np.mean(shorteset_path_lengths), color='red', label='True mean')
    # plt.title('Mean distance to allosteric pathway')
    # plt.xlabel('Shortest path length (sum of edge weights)')
    # plt.legend()
    # plt.savefig(f'results/allosteric_network_distance/mean_dist_to_pathway_{path_id}_{threshold_str}.png')

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Generate and save a graph from protein interactions')
    
    # Add command-line arguments
    parser.add_argument('-g', '--graph_path', type=str, help='path to the graph pickle file')
    parser.add_argument('-p', '--shortest_path_path', type=str, help='path to the shortest path file')
    parser.add_argument('-t', '--threshold', type=str, help='Min frequency to consider an edge in the shortest_path')
    parser.add_argument('-r', '--replicate', type=str, help='Replicate id to analyze')
    
    args = parser.parse_args()
    
    main(**vars(args))