DESC='''
MSC network properties
'''

import argparse
from itertools import product
import networkx as nx
import numpy as np
from os import environ,getenv
from os.path import basename
import pandas as pd
from pdb import set_trace
import scipy
import sqlite3
import sys

sys.path.append('../../../SPREAD_multipathway_simulator/simulator/scripts')
import msc_network as msc

pd.options.mode.chained_assignment = None

def main():
    # parser
    parser=argparse.ArgumentParser(description=DESC, 
            formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-i", "--network", required = True, help="MSC network")
    parser.add_argument("-r", "--range", required = True, type = float, help="range")
    parser.add_argument("-o", "--output_prefix", required = True, help="out prefix")
    args = parser.parse_args()

    # read network using msc
    network = msc.MultiScaleNet()
    network.read_from_file(args.network)
    network.display_summary()
    network.nodes[0] = network.nodes[0].set_index('node')
    network.edges[0] = network.edges[0][network.edges[0].moore <= args.range]
    hierarchy = network.hierarchy[network.hierarchy.parent != -1]

    # month
    sd_edges = network.edges[0][['source', 'target']]
    intra_loc_edges = hierarchy[
            ['parent', 'child']].groupby('parent').apply(
                    lambda x: nx.to_pandas_edgelist(
                        nx.complete_graph(x.child.tolist()))).reset_index()[
                                ['source', 'target']]

    with open(f'{args.output_prefix}_{args.range}.csv', 'w') as f:
        f.write('network,range,month,spectral_radius,spectral_radius_unweighted,diameter\n')
        for m in range(1,13):
            # short distance
            cells = network.nodes[0][f'{m}']
            sd_edges['weight'] = sd_edges.source.map(cells)

            # intra-locality
            intra_loc_edges['weight'] = intra_loc_edges.source.map(cells)

            # inter-locality
            Fld = network.edges[1][network.edges[1].month == m]
            meta_edge_list = []
            for ind, e in Fld.iterrows():
                x = hierarchy[hierarchy.parent == e.source].child.to_list()
                y = hierarchy[hierarchy.parent == e.target].child.to_list()
                meta_edge = pd.DataFrame.from_records(product(x,y)).rename(columns = {
                    0: 'source', 1:'target'})
                meta_edge['weight'] = meta_edge.source.map(cells) * e.weight
                meta_edge_list.append(meta_edge)
            inter_loc_edges = pd.concat(meta_edge_list)

            # combine networks
            edges = pd.concat([sd_edges, intra_loc_edges, inter_loc_edges]).groupby(
                ['source', 'target']).sum().reset_index()

            # diameter
            G = nx.from_pandas_edgelist(sd_edges, create_using = nx.DiGraph())
            diam = 0
            for comp in nx.strongly_connected_components(G):
                diam = max(diam, nx.diameter(G.subgraph(comp)))

            # eigen
            mat = scipy.sparse.csr_matrix((edges.weight, (edges.source, edges.target)))
            mat_unweighted = scipy.sparse.csr_matrix(
                    (np.ones(edges.shape[0]),(edges.source, edges.target)))
            eig = scipy.sparse.linalg.eigs(mat, k = 1)
            eig_unweighted = scipy.sparse.linalg.eigs(mat_unweighted, k = 1)
            
            f.write(f'{args.network},{args.range},{m},{eig[0][0].real},{eig_unweighted[0][0].real},{diam}\n')

if __name__== "__main__":
    main()
