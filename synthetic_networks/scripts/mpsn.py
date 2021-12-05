#!/usr/bin/python
DESC = '''Synthetic multipathway graph model

By: NP

Pending (AA):
    visualize: https://stackoverflow.com/questions/35109590/how-to-graph-nodes-on-a-grid-in-networkx
'''

import argparse
import csv
import logging
import math
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from numpy import linalg as LA
import scipy #Scipy 1.7? I think is incomplatible with networkx need to downgrade to scipy 1.4.1? to get it to work properly
import os
import pandas as pd
from pdb import set_trace
from scipy.sparse.linalg import eigs
from scipy.sparse import csgraph
#from scipy import stats
# from pandas import *: AA: avoid such 'import *' statements
import sqlite3
from sqlite3 import Error
import random
#import seaborn as sns #can't use problem with the stats from scipy so I can't use


FORMAT = "%(levelname)s:%(filename)s:%(funcName)s:%(message)s"


def self_mediated_dispersal(dispersalRange, side):
    H = nx.Graph()
    for i in range(0, side ** 2 - 1):
        col = i % side  # initializing which column its on
        row = i // side  # initializing rows
        # AA: check whether np.floor is correct
        for x in range(0, int(np.floor(dispersalRange)) + 1):  # is connecting the x direction columns
            for y in range(0, int(np.floor(dispersalRange)) + 1):  # is connecting the y direction rows
                rangeSquared = dispersalRange ** 2
                distanceSquared = x ** 2 + y ** 2
                if x == 0 and y == 0:
                    continue
                if col + x < side and row + y < side and distanceSquared <= rangeSquared:
                    H.add_edge((row, col), (row + y, col + x))
    return H


# AA: function names should be all lowercase (with underscores if required)
def locality_clique(gridSideLength, regionSideLength, localitySideLength, localityNumber):  # clique template
    # side is total number of localitySideLength
    # regionSideLength is number of localitySideLength in a regionSideLength
    # row is how many localitySideLength are in meta nodes 0 is the first meta node
    # regionSideLength number is which meta node it is/which regionSideLength going from left to right up to down
    x = (
                    regionSideLength - localitySideLength) // 2  # of meta node localitySideLength - number of white row nodes divided by 2. Gives the number of non white nodes surrounding white nodes in meta node
    # numberOfSquares = side**2//regionSideLength//regionSideLength
    initialRow = (localityNumber // (gridSideLength // regionSideLength)) * regionSideLength+x #Error where I forgot to include the x
    initialCol = (localityNumber % (gridSideLength // regionSideLength)) * regionSideLength+x
    G = nx.Graph()
    for i in range(0, localitySideLength):
        for j in range(0, localitySideLength):
            for k in range(0, localitySideLength):
                for l in range(0, localitySideLength):
                    if not (i == k and j == l):
                        G.add_edge((initialRow + i, initialCol + j), (initialRow + k,
                                                                      initialCol + l))  # iterates through every combonation of nodes and pairs them all togehter
                        G.add_edge((initialRow + k, initialCol + l), (initialRow + i,
                                                                      initialCol + j))  # probably don't need this line as in networkx edges are undirect

    return G  # returns a graph of all the white nodes in a metanode connected to every other white node


def locality_star(gridSideLength, regionSideLength, localitySideLength, localityNumber):
    midRow = localitySideLength // 2  # picks either the middle node or the bottom right of the smallest regionSideLength node for white localitySideLength
    midCol = localitySideLength // 2  # if even picks bottom right of 2x2 regionSideLength in the middle

    z = (
                    regionSideLength - localitySideLength) // 2  # of meta node localitySideLength - number of white row nodes divided by 2. Gives the number of non white nodes surrounding white nodes in meta node
    # numberOfSquares = gridSideLength*gridSideLength//regionSideLength//regionSideLength
    initialRow = (localityNumber // (gridSideLength // regionSideLength)) * regionSideLength + z
    initialCol = (localityNumber % (gridSideLength // regionSideLength)) * regionSideLength + z
    G = nx.Graph()
    for x in range(0, localitySideLength):
        for y in range(0, localitySideLength):
            if not (x == midRow and y == midCol):
                G.add_edge((initialRow + x, initialCol + y), (
                initialRow + midRow, initialCol + midCol))  # connects the head node to every other node except itself


    return G


def complete_bipartite(graph1,
                       graph2):  # function takes two graphs and returns a complete bipartate graph of one node set to another
    nodes1 = list(graph1.nodes)
    nodes2 = list(graph2.nodes)
    newGraph = nx.Graph()
    for x in range(0, len(nodes1)):
        for y in range(0, len(nodes2)):
            newGraph.add_edge(nodes1[x], nodes2[y])
    return newGraph

def heat_map(G,betweenness,filename,title):
    plt.clf()
    numNodes = len(list(G.nodes()))
    numRows = math.sqrt(numNodes)
    numRows = int(numRows)
    arr = [[0 for i in range(0,numRows)] for j in range(0,numRows)]
    for x in range(0,numRows):
        for y in range(0,numRows):
            arr[x][y]=betweenness[x,y]

    _ = plt.imshow(arr, cmap='autumn_r', interpolation='nearest')
    plt.title(title)
    plt.colorbar()
    plt.axis('off')
    #ax = sns.heatmap(arr, linewidth=0.5, cmap='coolwarm')
    plt.savefig(filename, bbox_inches=0)


def creating_files(G,GS,GL,GLD,localityNum,localityNodes,longDistanceEdges):
    #Need to check if simulation automatically assume bi-directionallity or not I don't think it does
    THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
    my_file = os.path.join(THIS_FOLDER, '0.nodes')
    current_name = os.path.join(THIS_FOLDER, 'nodes0.csv')
    # absolute path
    database = my_file
    nodesFile = '0.nodes'
    if os.path.exists(nodesFile):
        os.remove('0.nodes')

    # creating maps
    nodes = G.nodes()
    nodes = sorted(nodes, key=lambda tup: (tup[0], tup[1]))
    numRows = int(np.sqrt(len(nodes)))

    with open('0.nodes', 'w') as storage_file:
        storage_writer = csv.writer(storage_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
        data = ['node', 'node_x','node_y',1,2,3,4,5,6,7,8,9,10,11,12]
        storage_writer.writerow(data)
        for node in nodes:
            data = [node[0]*numRows+node[1],node[0],node[1],1,1,1,1,1,1,1,1,1,1,1,1]
            storage_writer.writerow(data)
        storage_file.close()
        #os.rename(current_name,my_file)

    edgeFile = '0.edges'
    if os.path.exists(edgeFile):
        os.remove('0.edges')
    with open('0.edges', 'w') as storage_file:
        storage_writer = csv.writer(storage_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
        data = ['source', 'target','moore','source_x','source_y', 'target_x','target_y']
        storage_writer.writerow(data)
        #print(G.edges())
        edges = list(G.edges())

        for edge in edges:
            node00 = edge[0][0]
            node01 = edge[0][1]
            node10 = edge[1][0]
            node11 = edge[1][1]
            nodeNumber0 = node00*numRows+node01
            nodeNumber1 = node10*numRows+node11
            distanceX = (node10-node00)**2
            distanceY = (node11-node01)**2
            distanceTotal = math.sqrt(distanceX+distanceY)
            data = [nodeNumber0,nodeNumber1,distanceTotal,edge[0][0],edge[0][1],edge[1][0],edge[1][1]]
            storage_writer.writerow(data)
            data = [nodeNumber1,nodeNumber0,distanceTotal,edge[1][0],edge[1][1],edge[0][0],edge[0][1]]
            storage_writer.writerow(data)
        storage_file.close()

    numberOfNodes = len(list(G.nodes()))
    nodesFile = '1.nodes'
    if os.path.exists(nodesFile):
        os.remove('1.nodes')
    with open('1.nodes', 'w') as storage_file:
        storage_writer = csv.writer(storage_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL,
                                    lineterminator='\n')
        data = ['node', 'locality_number']
        storage_writer.writerow(data)
        nodes = G.nodes()
        nodes = sorted(nodes, key=lambda tup: (tup[0], tup[1]))
        x = 0
        for node in range(0,localityNum):
            localityNumber = node+numberOfNodes+1000
            data = [localityNumber, node]
            storage_writer.writerow(data)
        storage_file.close()

    nodesFile='1.edges'
    if os.path.exists(nodesFile):
        os.remove('1.edges')
    with open('1.edges', 'w') as storage_file:
        storage_writer = csv.writer(storage_file, 
                delimiter=',', 
                quotechar='"', 
                quoting=csv.QUOTE_MINIMAL, 
                lineterminator='\n')
        data = ["source", "target", "weight", "month", "cell_id_x", "cell_id_y", "source_region", "target_region"]
        storage_writer.writerow(data)
        nodes = G.nodes()
        nodes = sorted(nodes, key=lambda tup: (tup[0], tup[1]))
        for edge in longDistanceEdges.edges():
            edge0 = edge[0]+numberOfNodes+1000
            edge1 = edge[1]+numberOfNodes+1000
            for y in range(1,13): #makes sure that it goes through every single month and that its constant throughout
                data = [edge0,edge1,1,y,0,0,edge[0],edge[1]]
                storage_writer.writerow(data)
                data = [edge1, edge0, 1, y, 0, 0, edge[1], edge[0]]
                storage_writer.writerow(data)
        storage_file.close()

    nodesFile = 'hierarchy.tree'
    if os.path.exists(nodesFile):
        os.remove('hierarchy.tree')

    seed_loc = pd.DataFrame(columns = ['node', 'probability', 'locality'])
    with open('hierarchy.tree', 'w') as storage_file:
        storage_writer = csv.writer(storage_file, 
                delimiter=',', 
                quotechar='"', 
                quoting=csv.QUOTE_MINIMAL, 
                lineterminator='\n')
        data = ['parent', 'child', 'parent_node', 'child_node_x', 'child_node_y']
        storage_writer.writerow(data)
        nodes = G.nodes()
        nodes = sorted(nodes, key=lambda tup: (tup[0], tup[1]))
        x = 0
        for node in range(0, localityNum):
            localityNumber = node + numberOfNodes + 1000
            data = [-1, localityNumber, -1, localityNum]
            storage_writer.writerow(data)
        for locality in range(0, localityNum):
            localityNumber = locality + numberOfNodes + 1000
            for node in localityNodes[locality]:
                nodeNumber = node[0] * numRows + node[1]
                data = [localityNumber, nodeNumber, locality, node[0], node[1]]
                seed_loc = seed_loc.append({
                    'node': nodeNumber,
                    'locality': locality}, ignore_index = True)
                storage_writer.writerow(data)
        onlyShortDistance = [node for node in G.nodes() if node not in GL.nodes()]
        #print(onlyShortDistance)
        for node in onlyShortDistance:
            nodeNumber = node[0] * numRows + node[1]
            data = [-1,nodeNumber,-1,nodes[0],node[1]]
            seed_loc = seed_loc.append({
                'node': nodeNumber, 
                'locality': -1}, ignore_index = True)
            storage_writer.writerow(data)
        storage_file.close()

    # seeding
    percentageSeedNodes = 5
    df = pd.DataFrame({'node': np.arange(G.number_of_nodes()),
        'probability': np.full(G.number_of_nodes(), percentageSeedNodes/100)})
    df.to_csv('seed_all.csv', index = False)
    numNodesInLocality = (seed_loc.locality != -1).sum()
    localityProb = min(1, G.number_of_nodes()/numNodesInLocality \
            * percentageSeedNodes / 100)
    nonLocalityProb = max(0, (G.number_of_nodes() * percentageSeedNodes / 100 - \
            numNodesInLocality * localityProb) / G.number_of_nodes())
    seed_loc.probability = nonLocalityProb
    seed_loc.loc[seed_loc.locality != -1, 'probability'] = localityProb
    seed_loc = seed_loc.astype({'node': int})
    seed_loc.to_csv('seed_loc.csv', index = False)

    ## df = pd.DataFrame({'node': nodesInLocality,
    ##     'probability': np.full(numNodesInLocality, localityProb}))
    ## df = df.append(pd.DataFrame({'node': nodesInLocal


def main():
    parser = argparse.ArgumentParser(description=DESC,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--number_of_nodes", required=True,
                        help="Number of nodes in the square grid. Must be a square.",
                        type=int)
    parser.add_argument("--number_of_regions", required=True,
                        help="Number of square regions in the grid. Must be a square and be a factor of number_of_nodes.",
                        type=int)
    parser.add_argument("--range", required=True,
                        help="Distance parameter", type=float)
    parser.add_argument("--locality_size", required=True,
                        help="Number of nodes in a locality, which forms a square grid within a region. Must be a square and have the same parity as number_of_nodes/number_of_regions.",
                        type=int)
    parser.add_argument("--locality_graph", required=True,
                        help="Type of locality graph (star/clique)")
    parser.add_argument("--long_distance_type", required=True,
                        help="Type of graph for long distance pathway (ER/CL/SF)")
    parser.add_argument("--ld_param", nargs="+",
                        help="Long Distance parameters",type=float)
    parser.add_argument("--seed", required=True, 
            help="Random seed", type=int)
    parser.add_argument("--directed", action="store_true",
                        help="says that the graph is a directional graph")
    parser.add_argument("-m", "--multi", action="store_true",
                        help="says that the graph is a multi-edged graph")
    parser.add_argument("--suppress_properties", action="store_true",
                        default = False,
                        help="Mode to just create graphs.")
    parser.add_argument("-d", "--debug", action="store_true")
    parser.add_argument("-q", "--quiet", action="store_true")

    args = parser.parse_args()

    # set logger
    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format=FORMAT)
    elif args.quiet:
        logging.basicConfig(level=logging.WARNING, format=FORMAT)
    else:
        logging.basicConfig(level=logging.INFO, format=FORMAT)

    # checking if constraints are satisfied
    if not np.sqrt(args.number_of_nodes).is_integer():
        raise ValueError('number of nodes must be a square')
    if not np.sqrt(args.number_of_regions).is_integer():
        raise ValueError('number of regions must be a square')
    if args.number_of_regions > args.number_of_nodes:
        raise ValueError('number of regions > number_of_nodes')
    if not np.sqrt(args.locality_size).is_integer():
        raise ValueError('locality size must be a square')
    if args.locality_size > args.number_of_nodes / args.number_of_regions:
        raise ValueError('locality size > number of nodes in a region')
    if (args.number_of_nodes / args.number_of_regions) % 2 != args.locality_size % 2:
        raise ValueError(f'locality size {args.locality_size} and number of regions {args.number_of_regions} must have the same parity')

    if args.seed:
        random.seed(args.seed)
        np.random.seed(args.seed)


    gridSideLength = int(np.sqrt(args.number_of_nodes))
    regionSideLength = int(np.sqrt(args.number_of_nodes // args.number_of_regions))
    localitySideLength = int(np.sqrt(args.locality_size))
    dispersalRange = args.range
    numberOfRegions = args.number_of_regions

    GS = self_mediated_dispersal(dispersalRange,
                                 gridSideLength)  # AA: no need to return G since you do not create a new object in this function

    localityGraphList = []  # holds individual graphs of localities
    # GLNodes = [0]*(gridSideLength**2//regionSideLength**2) #array of size number of meta nodes
    localityGraphsSet = []
    whichLocality = "Normal"  # Locality is set to square
    GL = nx.Graph()
    for locality in range(0, args.number_of_regions):
        if args.locality_graph == 'star':
            localityGraph = locality_star(gridSideLength,
                                          regionSideLength,
                                          localitySideLength, locality)  # Generates a star graph for locality
        else:             #elif args.locality_graph == 'clique':
            localityGraph = locality_clique(gridSideLength,
                                            regionSideLength,
                                            localitySideLength, locality)  # Generates a clique graph for locality
        GL = nx.compose(GL, localityGraph)
        localityGraphsSet.append(
            localityGraph)  # which nodes are the white nodes that will be connected via long distance
    # LocalityGraph = nx.compose(LocalityGraph, Glocality[locality]) #merge all the locality graphs together

    if args.long_distance_type == 'ER':
        prob = min(1,float(args.ld_param[0])/args.number_of_regions)  # get probability
        H = nx.erdos_renyi_graph(args.number_of_regions, prob)
    elif args.long_distance_type == 'CL':
        H = nx.expected_degree_graph(
            args.ld_param)  # generates the graph using the array of weights that are the input parameters
    elif args.long_distance_type == 'SF':
        H = nx.scale_free_graph(args.number_of_regions)  # generates a scale free graph

    GLD = nx.Graph()

    for edge in H.edges():
        if args.number_of_regions != 1 and edge[0] != edge[1]   :
            GLD = nx.compose(GLD,
                         complete_bipartite(localityGraphsSet[edge[0]],
                                            localityGraphsSet[edge[1]]))

    # Merge all graphs
    # AA: multigraph not working in many cases
    G = nx.Graph()
    G.add_edges_from(GS.edges, label="S")  # Merging natural spread
    G.add_edges_from(GL.edges, label="L")  # Merging Locality
    G.add_edges_from(GLD.edges, label="LD")  # Merging all the graphs together

    # compute measures
    # AA: CompleteGraph in graph theory means all edges are present
    ## adjMatrix = nx.to_numpy_array(G)  # convert to numpy array to find eigen values
    ## # AA: add 1st eigenvalue and 2nd eigenvalue
    ## adjSpectrum = LA.eigvals(adjMatrix)  # measures eigenvalue
    ## adjSpectrum.sort()
    adjMat = nx.to_scipy_sparse_matrix(G).asfptype()
    spectralRadius, eigvec = eigs(adjMat, k=1)  # measures eigenvalue
    spectralRadius = np.real(spectralRadius[0])
    #spectralRadius = adjSpectrum[len(adjSpectrum)-1]
    ## eigenValue2 = adjSpectrum[len(adjSpectrum)-2]
    # AA: also compute laplacian spectrum https://networkx.org/documentation/stable/reference/generated/networkx.linalg.laplacianmatrix.laplacian_matrix.html: again 1st and 2nd
    lapMat = csgraph.laplacian(adjMat, normed=False)
    laplacian, vec = eigs(lapMat, k=2, which='SM')
    laplacian = np.real(laplacian[1])

    creating_files(G,GS,GL,GLD,numberOfRegions,localityGraphsSet,H)

    if not args.suppress_properties:
        diameter = nx.diameter(G)
        betweennessUnweighted = nx.betweenness_centrality(G, weight=None)  # Measures betweenness
        maxCoreNum = max(nx.core_number(G).values())  # int
        heat_map(G,betweennessUnweighted,
                'betweenness_centrality.pdf',
                'Betweenness centrality')
        rowInd = range(int(np.sqrt(args.number_of_nodes)))
        x = np.meshgrid(rowInd,rowInd)
        betweennessEig = dict(zip(zip(x[1].flatten(),x[0].flatten()), 
            eigvec.real.transpose()[0]))
        heat_map(G,betweennessEig,
                'betweenness_eig.pdf',
                'Eigenvector centrality')
    # add all the data from the measures and properties of the graph into the SQL file.
    # AA: just need the print statement
    if not args.suppress_properties:
        print(f"INSERT OR REPLACE INTO mpsn_props (\
number_of_nodes,\
number_of_regions,\
range,\
locality_size,\
locality_graph,\
long_distance_type,\
epsilon,\
seed,\
number_of_edges,\
GS_number_of_edges,\
GL_number_of_edges,\
GLD_number_of_edges,\
max_core_number,\
spectral_radius,\
laplacian,\
diameter\
) VALUES (\
{args.number_of_nodes},\
{args.number_of_regions},\
{args.range},\
{args.locality_size},\
'{args.locality_graph}',\
'{args.long_distance_type}',\
{args.ld_param[0]},\
{args.seed},\
{G.number_of_edges()},\
{GS.number_of_edges()},\
{GL.number_of_edges()},\
{GLD.number_of_edges()},\
{maxCoreNum},\
{spectralRadius},\
{laplacian},\
{diameter}\
);")

if __name__ == '__main__':
    main()
