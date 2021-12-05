DESC='''Synthetic multipathway graph model

By: NP

'''

import argparse
import csv
import logging
import math
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from numpy import linalg as LA
import os
import pandas as pd
from pdb import set_trace
#from pandas import *: AA: avoid such 'import *' statements
import sqlite3
from sqlite3 import Error

FORMAT="%(levelname)s:%(filename)s:%(funcName)s:%(message)s"

def self_mediated_dispersal(dispersalRange,side):

    H=nx.Graph()
    for i in range(0,side**2-1):
        col = i % side #initializing which column its on
        row = i // side #initializing rows
        # AA: check whether np.floor is correct
        for x in range(0,int(np.floor(dispersalRange))+1): #is connecting the x direction columns
            for y in range (0,int(np.floor(dispersalRange))+1): #is connecting the y direction rows
                rangeSquared=dispersalRange**2
                distanceSquared = x**2+y**2
                if x==0 and y==0:
                    continue
                # AA: is x<=moore necessary?
                # AA: is row-y or col-x necessary?
                ## if col + x < side and row + y < side and x <= moore and y <= moore and not (x==0 and y == 0) and distance <= rangeSquared:
                ##     a.add_edge((row,col),(row+y,(col+x)))
                ## if col + x < side and row - y >= 0 and x <= moore and y <= moore and not (x==0 and y == 0) and distance <= rangeSquared:
                ##     a.add_edge((row,col),(row-y,(col+x)))
                ## if col - x >= 0 and row + y < side and x <= moore and y <= moore and not (x==0 and y == 0) and distance <= rangeSquared:
                ##     a.add_edge((row,col),(row+y,(col-x)))
                ## if col - x >= 0 and row - y >= 0 and x <= moore and y <= moore and not (x==0 and y == 0) and distance <= rangeSquared:
                ##     a.add_edge((row,col),(row-y,(col-x)))
                if col + x < side and row + y < side and distanceSquared <= rangeSquared:
                    H.add_edge((row,col),(row+y,col+x))
    return H

# AA: function names should be all lowercase (with underscores if required)
def locality_clique(gridSideLength, regionSideLength, localitySideLength, localityNumber): # clique template
    #side is total number of localitySideLength
    #regionSideLength is number of localitySideLength in a regionSideLength
    #row is how many localitySideLength are in meta nodes 0 is the first meta node
    #regionSideLength number is which meta node it is/which regionSideLength going from left to right up to down
    x = (regionSideLength - localitySideLength)//2 #of meta node localitySideLength - number of white row nodes divided by 2. Gives the number of non white nodes surrounding white nodes in meta node
    #numberOfSquares = side**2//regionSideLength//regionSideLength
    initialRow = (localityNumber // (gridSideLength//regionSideLength)) * regionSideLength
    initialCol = (localityNumber % (gridSideLength//regionSideLength)) * regionSideLength
    G = nx.Graph()
    for i in range(0,localitySideLength):
        for j in range(0,localitySideLength):
            for k in range(0,localitySideLength):
                for l in range (0,localitySideLength):
                    if not (i == k and j == l) :
                        G.add_edge((initialRow+i,initialCol+j),(initialRow+k,initialCol+l)) #iterates through every combonation of nodes and pairs them all togehter
                        G.add_edge((initialRow+k,initialCol+l),(initialRow+i,initialCol+j)) #probably don't need this line as in networkx edges are undirect

    return G #returns a graph of all the white nodes in a metanode connected to every other white node

def locality_star(gridSideLength, regionSideLength, localitySideLength, localityNumber):
    midRow = localitySideLength//2 #picks either the middle node or the bottom right of the smallest regionSideLength node for white localitySideLength
    midCol = localitySideLength//2 #if even picks bottom right of 2x2 regionSideLength in the middle

    z = (regionSideLength - localitySideLength)//2 #of meta node localitySideLength - number of white row nodes divided by 2. Gives the number of non white nodes surrounding white nodes in meta node
    #numberOfSquares = gridSideLength*gridSideLength//regionSideLength//regionSideLength
    initialRow = (localityNumber // (gridSideLength//regionSideLength)) * regionSideLength + z
    initialCol = (localityNumber % (gridSideLength//regionSideLength)) * regionSideLength + z
    G = nx.Graph()
    for x in range(0,localitySideLength):
        for y in range(0,localitySideLength):
            if not (x == midRow and y == midCol):
                G.add_edge((initialRow+x,initialCol+y), (initialRow+midRow,initialCol+midCol)) #connects the head node to every other node except itself
    return G

def complete_bipartite(graph1, graph2):#function takes two graphs and returns a complete bipartate graph of one node set to another
    nodes1 = list(graph1.nodes)
    nodes2 = list(graph2.nodes)
    newGraph = nx.Graph()
    for x in range (0,len(nodes1)):
        for y in range(0,len(nodes2)):
            newGraph.add_edge(nodes1[x], nodes2[y])
    return newGraph

def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)

    return conn


def update_task(conn, task):
    """
    update priority, begin_date, and end date of a task
    :param conn:
    :param task:
    :return: project id
    """
    #sql variable is defunct
    sql = ''' UPDATE projects
              SET   localityType = ?,
                    diameter = ?,
                    mooreRange = ?,
                    rowLength = ?,
                    squareRowSize = ?,
                    whiteNodeRows = ?,
                    LDGraphType = ?,
                    maxCoreNumber = ?,
                    spectralRadius = ?'''
    #inserting into which columns in the sql file order of variables and columns.
    sqlite_insert_with_param = """INSERT INTO projects
                          (localityType, diameter, mooreRange, rowLength, squareRowSize, whiteNodeRows, LDGraphType, maxCoreNumber, spectralRadius)
                          VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);"""
    #inserting into sql file.
    cur = conn.cursor()
    cur.execute(sqlite_insert_with_param, task)
    conn.commit()

def create_table(conn, create_table_sql):
    """ create a table from the create_table_sql statement
    :param conn: Connection object
    :param create_table_sql: a CREATE TABLE statement
    :return:
    """
    try:
        #Creates a table if table is not already created.
        c = conn.cursor()
        c.execute(create_table_sql)
        conn.commit()
        c.close()
    except Error as e:
        print(e)


def main():
    parser=argparse.ArgumentParser(description=DESC,
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--number_of_nodes",required=True,
            help="Number of nodes in the square grid. Must be a square.",
            type=int)
    parser.add_argument("--number_of_regions",required=True,
            help="Number of square regions in the grid. Must be a square and be a factor of number_of_nodes.",
            type=int)
    parser.add_argument("--range", required=True,
            help="Distance parameter", type=float)
    parser.add_argument("--locality_size",required=True,
            help="Number of nodes in a locality, which forms a square grid within a region. Must be a square and have the same parity as number_of_nodes/number_of_regions.",
            type=int)
    parser.add_argument("--locality_graph",required=True,
            help="Type of locality graph (star/clique)")
    parser.add_argument("--long_distance_type", required=True,
            help="Type of graph for long distance pathway (ER/CL/SF)")
    parser.add_argument("--ld_param",nargs="+", 
            help="Long Distance parameters")
    parser.add_argument("--directed",action="store_true",
            help="says that the graph is a directional graph")
    parser.add_argument("-m","--multi",action="store_true", 
            help="says that the graph is a multi-edged graph")
    parser.add_argument("-d", "--debug", action="store_true")
    parser.add_argument("-q", "--quiet", action="store_true")

    args = parser.parse_args()

    # set logger
    if args.debug:
       logging.basicConfig(level=logging.DEBUG,format=FORMAT)
    elif args.quiet:
       logging.basicConfig(level=logging.WARNING,format=FORMAT)
    else:
       logging.basicConfig(level=logging.INFO,format=FORMAT)

    # checking if constraints are satisfied
    if not np.sqrt(args.number_of_nodes).is_integer():
        raise ValueError('number_of_nodes must be a square')
    if not np.sqrt(args.number_of_regions).is_integer():
        raise ValueError('number_of_regions must be a square')
    if args.number_of_regions > args.number_of_nodes:
        raise ValueError('number_of_regions > number_of_nodes')
    if not np.sqrt(args.locality_size).is_integer():
        raise ValueError('locality_size must be a square')
    if args.locality_size > args.number_of_nodes/args.number_of_regions:
        raise ValueError('locality_size > number of nodes in a region')
    if args.number_of_nodes/args.number_of_regions % 2 != args.locality_size % 2:
        raise ValueError('locality_size and region size must have the same parity')

    gridSideLength = int(np.sqrt(args.number_of_nodes))
    regionSideLength = int(np.sqrt(args.number_of_nodes//args.number_of_regions))
    localitySideLength = int(np.sqrt(args.locality_size))
    dispersalRange = args.range
    numberOfRegions=args.number_of_regions

    GS=self_mediated_dispersal(dispersalRange,gridSideLength)    # AA: no need to return G since you do not create a new object in this function

    localityGraphList = []  # holds individual graphs of localities
    #GLNodes = [0]*(gridSideLength**2//regionSideLength**2) #array of size number of meta nodes
    localityGraphsSet=[]
    whichLocality = "Normal" #Locality is set to square
    GL = nx.Graph()
    for locality in range(0,args.number_of_regions):
        if args.locality_graph=='star':
            localityGraph=locality_star(gridSideLength, 
                regionSideLength, 
                localitySideLength,locality) # Generates a star graph for locality
        elif args.locality_graph=='clique':
            localityGraph=locality_clique(gridSideLength, 
                regionSideLength, 
                localitySideLength,locality)  # Generates a clique graph for locality
        GL=nx.compose(GL,localityGraph)
        localityGraphsSet.append(localityGraph) #which nodes are the white nodes that will be connected via long distance
    #LocalityGraph = nx.compose(LocalityGraph, Glocality[locality]) #merge all the locality graphs together

    if args.long_distance_type=='ER':
        prob=float(args.ld_param[0])    # get probability
        H = nx.erdos_renyi_graph(args.number_of_regions,prob)
    elif args.long_distance_type=='CL':
        H = nx.expected_degree_graph(args.Chung) #generates the graph using the array of weights that are the input parameters
    elif args.long_distance_type=='SF': 
        H = nx.scale_free_graph(args.number_of_regions) #generates a scale free graph
    
    GLD = nx.Graph()
    for edge in H.edges():
        GLD = nx.compose(GLD,
            complete_bipartite(localityGraphsSet[edge[0]],
                localityGraphsSet[edge[1]]))
    
    # Merge all graphs
    # AA: multigraph not working in many cases
    G = nx.Graph()
    G.add_edges_from(GS.edges, label="S") #Merging natural spread
    G.add_edges_from(GL.edges, label="L") #Merging Locality
    G.add_edges_from(GLD.edges, label="LD") #Merging all the graphs together

    # compute measures
    # AA: CompleteGraph in graph theory means all edges are present
    adjMatrix = nx.to_numpy_array(G) #convert to numpy array to find eigen values
    # AA: add 1st eigenvalue and 2nd eigenvalue
    adjSpectrum = LA.eigvals(adjMatrix).sort() #measures eigenvalue
    # AA: also compute laplacian spectrum https://networkx.org/documentation/stable/reference/generated/networkx.linalg.laplacianmatrix.laplacian_matrix.html: again 1st and 2nd
    graphDiameter = nx.diameter(G) #Need to fix finding graph diameter apparently when I add the layers it has trouble figuring out that the graph is complete
    betweennessUnweighted = nx.betweenness_centrality(G,weight=None) #Measures betweenness

    # AA: the below two betweenness are used for weighted graphs.
    #betweennessInvWeighted = nx.betweenness_centrality(G,weight='invweight') #Measures betweensness but invweighted
    #betweennessNegLogWeighted = nx.betweenness_centrality(G,weight='invlogweight') #measures betweeness negative log weight
    maxCoreNum=max(nx.core_number(G).values()) #int
    np.set_printoptions(threshold=np.inf)

    #add all the data from the measures and properties of the graph into the SQL file.
    # AA: just need the print statement
    print(f'INSERT INTO projects (\
number_of_nodes,\
number_of_regions,\
range,\
locality_size,\
locality_graph,\
long_distance_type,\
max_core_number,\
spectral_radius,\
diameter\
) VALUES (\
{args.number_of_nodes},\
{args.number_of_regions},\
{args.range},\
{args.locality_size},\
{args.locality_graph},\
{args.long_distance_type},\
{maxCoreNum},\
{spectralRadius},\
{diameter}\
)')

if __name__=='__main__':
    main()
