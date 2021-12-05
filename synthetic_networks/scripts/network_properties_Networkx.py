#!/usr/bin/env python
###########################################################################
# AA 2016-04-26 - created
# AA 2016-05-11 - help & argument passing documented
# Purpose: Constructs csv files compiling data and network analysis given
# a set of input files.
###########################################################################

###########################################################################
# import
###########################################################################
import argparse
import logging
import networkx as nx
import pdb
import sys
import os
from operator import itemgetter
from math import log

########## constants
GLOBAL_TRAITS_FILE="global_traits.csv"
MAX_CORE_FILE="max_core.txt"

DESC="""
This script takes as input the following: 
1. the directed commodity flow/pathway network with nodes as locations,
2. (optional) threshold value to filter edges with weight less than threshold.
All output files will written to the specified/default output folder.
Many output files are of the form ./out/traits_[indicator].csv files for
each network/node trait studied. Each file is a table sorted in descending
order by the indicator in its name. In addition, it dumps the global traits
file (%s) and the set of nodes in the max core (%s).

""" %(GLOBAL_TRAITS_FILE,MAX_CORE_FILE)

## # calculates the strength of a given node, calculated as the total volume of
## # exports
## def strength(graph, node):
##    total = 0
##    neighbors = G[node]
##    for neigh in neighbors:
##       total += neighbors[neigh]['weight']
##    return total

###########################################################################
#Exports specific data sorted by one indicator, and outputs to a .csv 
#Inputs
#   file     - file object to store the output on
#   datalist - list of indicator names, which are written as column
#         heads (e.g. ["Strength", "Degree"])
#    sortOn     - Indicator which the output table will be sorted with
#         respect to.
#   data    - array of arrays of the data. data[i] is a list of all 
#         data relating to one trait mentioned in datalist.
###########################################################################
def export(file, datalist, sortOn, data):
   titleStr = ""
   with open(file, "w") as File:
      ## for title in datalist:
      ##    titleStr += title + ","
      ## titleStr = titleStr[:-1] + "\n"
      titleStr = "%s,%s\n" %(datalist[0],datalist[sortOn])
      File.write(titleStr)
      
      data = sorted(data, key = itemgetter(sortOn,0),reverse=True)
      for point in data:
      ##    dataStr = ",".join([str(i) for i in point])
      ##    dataStr += "\n"
         dataStr="%s,%s\n" %(str(point[0]),str(point[sortOn]))
         File.write(dataStr)

###########################################################################
# Generates the csv scripts as specified by command line arguments.
###########################################################################
if __name__=='__main__':
   parser = argparse.ArgumentParser(description=DESC,formatter_class=argparse.RawTextHelpFormatter)
   parser.add_argument("network",help="the network being analyzed")  
   parser.add_argument("-t","--threshold", action="store",type=float,help="any edge with weight below the threshold will be ignored",default=0)
   parser.add_argument("-o","--output_folder", action="store",help="output folder name",default="out")
   parser.add_argument("-u","--undirected", action="store_true",help="undirected")
   parser.add_argument("-d","--delimiter", action="store",help="delimiter (default is whitespace",default=" ")
   parser.add_argument("-v", "--verbose", action="store_true", help="print status messages to stdout")
   args = parser.parse_args()

   # set up main logging
   if args.verbose:
      logging.basicConfig(format='%(levelname)s: %(name)s: %(funcName)s: %(message)s', filemode='w', level=logging.INFO)
   else:
      logging.basicConfig(format='%(levelname)s: %(message)s', filemode='w', level=logging.WARNING)

   # load network
   if args.undirected:
      H=nx.Graph()
   else:
      H=nx.DiGraph()

   # handling different ways in which weight is specified
   try:
      G = nx.read_edgelist(args.network, create_using=H, delimiter=args.delimiter)
   except:
      G = nx.read_edgelist(args.network, create_using=H, delimiter=args.delimiter, data=(('weight',int),))

   # filter edges below threshold
   for e in G.edges():
      if G.edge[e[0]][e[1]]['weight']<args.threshold:
         G.remove_edge(e[0],e[1])

   Gun = G.to_undirected()

   ###########################################################################
   # Assembles dictionaries for each trait being evaluated
   ###########################################################################
   nodes = G.nodes()
   degrees = G.degree(nodes)
   if not args.undirected:
      in_degrees = G.in_degree(nodes)
      out_degrees = G.out_degree(nodes)
      weighted_in_degrees = G.in_degree(nodes,weight='weight')
      weighted_out_degrees = G.out_degree(nodes,weight='weight')
   else:
      in_degrees = {}
      out_degrees = {}
      weighted_in_degrees = {}
      weighted_out_degrees = {}

   # Assigns new attributes to each node
   for node in nodes:
      for neigh,attr in G[node].items():
         attr['invweight'] = 1/attr['weight']
         attr['invlogweight'] = -log(attr['weight'])
   
   ## # generate dictionary of node strengths
   ## strengths = dict()
   ## for n in nodes:
   ##    strengths[n] = strength(G,n)

   # Calculates betweenness centralities
   betweennessUnweighted = nx.betweenness_centrality(G,weight=None)
   betweennessInvWeighted = nx.betweenness_centrality(G,weight='invweight')
   betweennessNegLogWeighted = nx.betweenness_centrality(G,weight='invlogweight')
   if not args.undirected:
      Gun.remove_edges_from(Gun.selfloop_edges())
   else:
      Gun=G
   CoreNum=nx.core_number(Gun)
   maxCoreNum=max(CoreNum.values())
   
   maxCore=[];
   for node in CoreNum:
      if CoreNum[node]==maxCoreNum:
         maxCore.append(node)

   # Add new traits to this   
   _traits = {"node":         nodes, 
      "degree":         degrees, 
      "bet_cent_inv_wt":   betweennessInvWeighted,
      "bet_cent_neg_log_wt":   betweennessNegLogWeighted,
      "bet_cent_no_wt":   betweennessUnweighted,
   }
   
   if not args.undirected:
      _traits["in_degree"]=in_degrees
      _traits["out_degree"]=out_degrees 
      _traits["weighted_in_degree"]=weighted_in_degrees
      _traits["weighted_out_degree"]=weighted_out_degrees

   titles = _traits.keys()
   traits = _traits.values()
   
   ###########################################################################
   # Creates and exports data tables
   ###########################################################################
   # Generates the data tables for exporting.
   data = []
   for i,n in enumerate(nodes):
      #Include any additional traits in this list
      row = [n.replace(",","")]
      
      for trait in traits[1:]:
         row.append(trait[n])
      data.append(row)

   # create output folder
   outputFolder=args.output_folder
   if not os.path.exists(outputFolder):
      os.makedirs(outputFolder)
   logging.info('Output folder: %s' %outputFolder)

   # Generates the node traits csv file
   for i in range(1,len(titles)):
      logging.info('Writing the %s file ...' %titles[i])
      export(outputFolder+'/traits_' + titles[i] + '.csv', titles, i, data)

   # Global traits
   logging.info('Writing global traits file ...')
   with open(outputFolder+'/'+GLOBAL_TRAITS_FILE,'w') as f:
      f.write('Number of nodes, %d\n' %G.number_of_nodes())
      f.write('Number of edges, %d\n' %G.number_of_edges())
      f.write('Max. core number, %d\n' %maxCoreNum)
      if args.undirected:
         f.write('Components, %d\n' %nx.number_connected_components(G))
         f.write('Diameter, %d\n' %nx.diameter(G))
         ecc = nx.eccentricity(G, sp=nx.all_pairs_dijkstra_path_length(G))
         f.write('Weighted diameter, %d\n' %nx.diameter(G,ecc))

   # core
   with open(outputFolder+'/'+MAX_CORE_FILE,'w') as f:
      for node in maxCore:
         f.write('%s\n' %str(node))
   logging.info('Done ...')

