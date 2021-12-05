import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy import linalg as LA
from pandas import *
import warnings
import argparse

print("hello")
sideLength = 15 #K1 #of rows and columns
Mrange = 2
SquareSide = 5 #K2 #of rows and columns in the meta nodes. K1 is always a multipleof K2
whiteRows = 3 #K3 #of rows for the white notes in the cneter of the meta nodes K3 < K2 and K3 % 2 == K2 % 2
parser = argparse.ArgumentParser() #initializing argparse
parser.add_argument("rowNumber", help="number of rows in grid", type=int) #Declaring which arguments I want
parser.add_argument("MooreRange", help="Moore Range", type=int)
parser.add_argument("SquareSize", help="number of rows in the meta node", type=int)
parser.add_argument("whiteNodeSize", help="number of rows in white nodes", type=int)
args = parser.parse_args()
sideLength = args.rowNumber #Setting arguments taken in when I run the file
Mrange = args.MooreRange
SquareSide = args.SquareSize
whiteRows = args.whiteNodeSize
#print(sideLength)
#print(Mrange)
#print(SquareSide)
#print(whiteRows)

def mooreRange(a, moore, side): #old function not in use for array based implementation

    for i in range(0,side*side-1):
        col = i % side #initializing which column its on
        row = i / side #initializing rows
        for x in range(0,moore+1): #is connecting the y direction
            for y in range (0, moore+1): #is connecting the x direction
                if col - y >= 0 and col + y < side and x + y <= moore and not (x == 0 and y == 0):
                    if row - x >= 0 and row+x<side:
                        a[i][i+y+side*x] = 1
                        a[i][i-y+side*x] = 1
                        a[i][i-y-side*x] = 1
                        a[i][i+y-side*x] = 1
                        a[i+y+side*x][i] = 1
                        a[i-y+side*x][i] = 1
                        a[i-y-side*x][i] = 1
                        a[i+y-side*x][i] = 1
                if col - y >= 0 and row + x < side and not (x == 0 and y == 0) and x + y <= moore:
                    a[i][i-y+x*side] = 1
                    a[i-y+x*side][i] = 1
                if col +y < side and row + x < side and not (x == 0 and y == 0) and x + y <= moore:
                    a[i][i+y+x*side] = 1
                    a[i+y+x*side][i] = 1
                if col - y >= 0 and row - x >= 0 and not (x == 0 and y == 0) and x + y <= moore:
                    a[i][i-y-x*side] = 1
                    a[i-y-x*side][i] = 1
                if col + y < side and row - x >= 0  and not (x == 0 and y == 0) and x + y <= moore:
                    a[i][i+y-x*side] = 1
                    a[i+y-x*side][i] = 1

    return a


def mooreRangeGraph(a, moore, side): #connects all the nodes in a grid pattern with the option for adjusting the moore range
    #a is the empty graph (probably don't need it)
    #moore is what the moore range is
    #side is how many rows/columns there are
    for i in range(0,side*side-1):
        col = i % side #initializing which column its on
        row = i // side #initializing rows
        for x in range(0,moore+1): #is connecting the x direction columns
            for y in range (0, moore+1): #is connecting the y direction rows
                """
                Changing all of this because it didn't work so I rewrote it. It didn't work cause I made a stupid mistake probably. New code much better though\


                if col - y >= 0 and col + y < side and x + y <= moore and not (x == 0 and y == 0): #Checking boundary conditions Make sure don't connect node to itself and doesn't go out of bounds
                    #distance from original node always has to be less than the moore range
                    if row - x >= 0 and row+x<side: #more checking boundary conditions
                        #I do integer division for row and modulus for column. I do this on the math that I do to get which node is which. Calculate it based on 0-number of nodes scale.
                        a.add_edge((i//12,i%12),((i+y+side*x)//12, (i+y+side*x)%12))
                        a.add_edge((i//12,i%12),((i-y+side*x)//12, (i-y+side*x)%12))
                        a.add_edge((i//12,i%12),((i-y-side*x)//12, (i-y-side*x)%12))
                        a.add_edge((i//12,i%12),((i-y-side*x)//12, (i-y-side*x)%12))
                        a.add_edge(((i+y+side*x)//12, (i+y+side*x)%12),(i//12,i%12))
                        a.add_edge(((i-y+side*x)//12, (i-y+side*x)%12),(i//12,i%12))
                        a.add_edge(((i-y-side*x)//12, (i-y-side*x)%12),(i//12,i%12))
                        a.add_edge(((i+y-side*x)//12, (i+y-side*x)%12),(i//12,i%12))


                if col - y >= 0 and row + x < side and not (x == 0 and y == 0) and x + y <= moore:

                    a.add_edge((i//12,i%12),((i-y+x*side)//12, (i-y+x*side)%12))
                    a.add_edge(((i-y+side*x)//12, (i-y+side*x)%12),(i//12,i%12))
                if col +y < side and row + x < side and not (x == 0 and y == 0) and x + y <= moore:

                    a.add_edge((i//12,i%12), ((i+y+side*x)//12, (i+y+x*side)%12))
                    a.add_edge(((i+y+side*x)//12,(i+y+side*x)%12),(i//12,i%12))
                if col - y >= 0 and row - x >= 0 and not (x == 0 and y == 0) and x + y <= moore:

                    a.add_edge((i//12,i%12),((i-y-side*x)//12, (i-y-x*side)%12))
                    a.add_edge(((i-y-side*x)//12, (i-y-side*x)%12),(i//12,i%12))
                if col + y < side and row - x >= 0  and not (x == 0 and y == 0) and x + y <= moore:

                    a.add_edge((i//12,i%12),((i+y-side*x)//12, (i+y-x*side)%12))
                    a.add_edge(((i+y-side*x)//12,(i+y-side*x)%12),(i//12,i%12))
                """

                #this checks is each of the four directions where it could go out of bounds and if its not going out of bounds it adds and edge to the graph
                if col + x < side and row + y < side and x + y <= moore and not (x==0 and y == 0):
                    a.add_edge((row,col),(row+y,(col+x)))
                if col + x < side and row - y >= 0 and x + y <= moore and not (x==0 and y == 0):
                    a.add_edge((row,col),(row-y,(col+x)))
                if col - x >= 0 and row + y < side and x + y <= moore and not (x==0 and y == 0):
                    a.add_edge((row,col),(row+y,(col-x)))
                if col - x >= 0 and row - y >= 0 and x + y <= moore and not (x==0 and y == 0):
                    a.add_edge((row,col),(row-y,(col-x)))


    return a




def connecting(L1, L2, a): #old function not in use
    for i in range (0,4):
        for j in range (0,4):
            for x in range (0,4):
                for y in range (0,4):
                    a[L1+i*12+j][L2+x*12+y] = 1
                    a[L2+x*12+y][L1+i*12+j] = 1

    return a
def connect(L1, L2, graph): #old function not using right now. Connects all the meta nodes together
    for i in range (0,4):
        for j in range (0,4):
            for x in range (0,4):
                for y in range (0,4):
                    graph.add_edge(L1+i*12+j, L2+x*12+y)
                    graph.add_edge(L2+x*12+y, L1+i*12+j)

    return graph

def locality(side, square, rows, squareNumber): #function returns every node in a locality connected to every other node in the locality
    #Side is row number
    #Square is number of squares
    #row is how many rows are in meta nodes 0 is the first meta node
    #square number is which meta node it is
    x = (square - rows)//2 #of meta node rows - number of white row nodes divided by 2. Gives the number of non white nodes surrounding white nodes in meta node
    numberOfSquares = side*side//square//square
    initialRow = (squareNumber // (side//square)) * square + x
    initialCol = (squareNumber % (side//square)) * square + x
    g = nx.Graph()
    for i in range(0,rows):
        for j in range(0,rows):
            for k in range(0,rows):
                for l in range (0,rows):
                    if not (i == k and j == l) :
                        g.add_edge((initialRow+i,initialCol+j),(initialRow+k,initialCol+l)) #iterates through every combonation of nodes and pairs them all togehter
                        g.add_edge((initialRow+k,initialCol+l),(initialRow+i,initialCol+j)) #probably don't need this line as in networkx edges are undirect

    return g #returns a graph of all the white nodes in a metanode connected to every other white node
def localityNodes(side, square, rows, squareNumber): #function that returns the nodes in the specific meta node region.
    #Side is row number
    #Square is number of squares
    #row is how many rows are in meta nodes 0 is the first meta node
    #square number is which meta node it is
    x = (square - rows)//2 #of meta node rows - number of white row nodes divided by 2. Gives the number of non white nodes surrounding white nodes in meta node
    numberOfSquares = side*side//square//square
    initialRow = (squareNumber // (side//square)) * square + x #Initial row for first node in meta node
    initialCol = (squareNumber % (side//square)) * square + x  #initial column for first node in meta node
    g = nx.Graph()
    for i in range(0,rows):
        for j in range(0,rows):
            g.add_node((initialRow+i,initialCol+j)) #adds each node in meta node to graph g


    return g #returns graph g

def completeBipartate(graph1, graph2):#function takes two graphs and returns a complete bipartate graph of one node set to another
    nodes1 = list(graph1.nodes)
    nodes2 = list(graph2.nodes)
    newGraph = nx.Graph()
    for x in range (0,len(nodes1)):
        for y in range(0,len(nodes2)):
            newGraph.add_edge(nodes1[x], nodes2[y])
    return newGraph




warnings.filterwarnings("ignore", category=UserWarning)


#OLD CODE NOT IN USE
#G = nx.Graph()
#my_rows, my_cols = (12,12)
#arr = [[0]*(sideLength*sideLength) for _ in range(sideLength*sideLength)]

#H = nx.path_graph(my_rows*my_cols)
#G.add_nodes_from(H)
#OLD CODE FROM SETTING UP ADJACENCY MATRIX NOT CURRENTLY IN USE
#i = 0
#j = 0
#while i < sideLength-1:
    #j = 0
    #while j < sideLength:

        #G.add_edge(i*sideLength+j,i*sideLength+j+1)

        #print(i*my_rows+j)
        #print(i*my_rows+j+1)
        #arr[i*my_rows+j][i*my_rows+j+1] = 1
        #arr[i*my_rows+j+1][i*my_rows+j] = 1



        #G.add_edge(i*sideLength+j,i*sideLength+sideLength+j)
        #arr[i*my_rows+j][i*my_rows+my_cols+j] = 1
        #arr[i*my_rows+my_cols+j][i*my_rows+j] = 1
        #j+=1
    #i+=1

#i = 0
#while i < sideLength*sideLength-1:
    #arr[i][i+1] = 1
    #arr[i+1][i] = 1
    #if i < 144-12:
        #arr[i][i+12] = 1
        #arr[i+12][i] = 1
    #i+=1


#arr = connecting(0,52,arr)
#arr = connecting(4,52,arr)
#arr = connecting(8,48,arr)
#arr = connecting(48,96,arr)
#arr = connecting(52,96,arr)
#arr = connecting(52,104,arr)
#arr = connecting(56,100,arr)
#arr = connecting(96,100,arr)
#arr = mooreRange(arr, Mrange, sideLength)
#G = connect(0,52,G)
#G = connect(4,52,G)
#G = connect(8,48,G)
#G = connect(48,96,G)
#G = connect(52,96,G)
#G = connect(52,104,G)
#G = connect(56,100,G)
#G = connect(96,100,G)

#nx.draw_random(G, with_labels = True)
#plt.savefig("Graph.png")
#pd.set_option('display.max_rows', None)
#pd.set_option('display.max_columns', None)
#pd.set_option('display.width', None)
#pd.set_option('display.max_colwidth', None)
#print(list(G.edges))
#print(DataFrame(arr))
#np.set_printoptions(threshold=np.inf)
#print(np.array(arr))
#eigen = LA.eigvals(arr)
#minEigen = eigen
#print(minEigen)
#a = arr
#minArray = arr
#index1 = 0
#index2 = 0
#for i in range (0,144):
    #for j in range (0,144):
        #if a[i][j] != 0:
            #old = a[i][j]
            #a[i][j] = 0
            #a[j][i] = 0
            #eigen = LA.eigvals(arr)
            #print(eigen[0])
            #a[i][j] = old
            #a[j][i] = old
        #if  eigen[0] < minEigen[0]:
            #minEigen = eigen
            #minArray = a
            #index1 = i
            #index2 = j
        #a = arr
    #print(i)
#print(minEigen)
#print(minArray)
#print(index1)
#print(index2)
#D = nx.diameter(G)
#print(D)
#print(arr[1][13])
#print(arr[2][14])



#G = nx.grid_2d_graph(sideLength, sideLength)

#Start of new good code
G = nx.Graph()
G = mooreRangeGraph(G,Mrange, sideLength) #initializing graph of natural spread with the possibility to adjust the moore range
#print(list(G.edges))
nx.draw_random(G, with_labels = True) #Plotting graph
plt.savefig("Graph.png") #plotting graph

Glocality = [0]*(sideLength*sideLength//SquareSide//SquareSide) #array of size number of meta nodes
GLNodes = [0]*(sideLength*sideLength//SquareSide//SquareSide) #array of size number of meta nodes

LocalityGraph = nx.Graph()
for x in range(0,sideLength*sideLength//SquareSide//SquareSide):
    Glocality[x] = locality(sideLength, SquareSide, whiteRows, x) #connecting all of the locality nodes to one another
    GLNodes[x] = localityNodes(sideLength, SquareSide, whiteRows, x) #which nodes are the white nodes that will be connected via long distance
    LocalityGraph = nx.compose(LocalityGraph, Glocality[x]) #merge all the locality graphs together


#Making all the long distance nodes into 1 graph
P = completeBipartate(GLNodes[0], GLNodes[4]) #creates complete Bipartate graph of the two node sets
#print(list(P.edges))
#print("test")
P = nx.compose(P, completeBipartate(GLNodes[1], GLNodes[4])) #Merging the new long distance with the old graph. Not necessarily a complete graph
P = nx.compose(P, completeBipartate(GLNodes[2], GLNodes[3]))
P = nx.compose(P, completeBipartate(GLNodes[3], GLNodes[6]))
P = nx.compose(P, completeBipartate(GLNodes[4], GLNodes[6]))
P = nx.compose(P, completeBipartate(GLNodes[4], GLNodes[8]))
P = nx.compose(P, completeBipartate(GLNodes[5], GLNodes[7]))
P = nx.compose(P, completeBipartate(GLNodes[6], GLNodes[7]))

#print(list(GLNodes[0].nodes))
#P = nx.cartesian_product(Glocality[0], Glocality[1])
#print(len(list(P.edges)))

CompleteGraph = nx.MultiGraph()
CompleteGraph.add_edges_from(P.edges, label="LongDistance") #Merging all the graphs together
#print(list(CompleteGraph.edges))
CompleteGraph.add_edges_from(LocalityGraph.edges, label="Locality") #Merging Locality
#print(list(CompleteGraph.edges))
CompleteGraph.add_edges_from(G, label="Natural") #Merging natural spread
#print(list(CompleteGraph.edges))
#print(list(LocalityGraph.edges))
#print(nx.get_edge_attributes(CompleteGraph,'label'))
nx.draw_random(G, with_labels = True) #Plotting graph
plt.savefig("Graph2.png") #plotting graph

#adjMatrix = nx.adjacency_matrix(CompleteGraph)
adjMatrix = nx.to_numpy_array(CompleteGraph) #convert to numpy array to find eigen values
eigenValue = LA.eigvals(adjMatrix)
#graphDiameter = nx.diameter(CompleteGraph) #Need to fix finding graph diameter apparently when I add the layers it has trouble figuring out that the graph is complete
print(eigenValue)
#print(graphDiameter)

