#!/usr/bin/python

from __future__ import print_function
from accessibility_lib import *
import sys,argparse


#Reading and treating arguments ----------------------------------

import argparse
parser = argparse.ArgumentParser(description="Calculate the accessibility value of each node of a network given as input.")
parser.add_argument("inputfile", help="Network file in .gml format as input.")
parser.add_argument("-o","--outputfile", help="File to output the given network with Accessibility values as a vertex attribute. Default: print values computed on stdout.",type=str)
parser.add_argument("-u","--undirected", help="Consider network as undirected.",action="store_true")
parser.add_argument("-w","--walk", help="Size of each random walk. Default=3",type=int,default=3)
parser.add_argument("-s","--steps", help="Number of random walks performed for each node. Default=200",type=int,default=200)
parser.add_argument("-n","--normalized", help="Make the accessibility values normalized.", action="store_true")
parser.add_argument("-v","--verbose", help="Make the software more verbose.",action="store_true")
args = parser.parse_args()



# ----------------------------------------------------------------

if args.verbose:
	print("Reading network from file \""+args.inputfile+"\"... ",end="")
g = Graph.Read_GML(args.inputfile)
if args.verbose:
	print("Done")

directed = g.is_directed()
if args.undirected:
	directed = False

N = len(g.vs)	
accessibility = np.zeros(shape=(N))

if args.normalized:
	for i in xrange(N):
		if args.verbose:
			print('Calculating progress: '+str(int(float(i+1)*100/N))+'%',end='\r')
			sys.stdout.flush()
		accessibility[i] = calculateNormalizedNodeDiversity(g,args.walk,i,args.steps,directed)
else:
	for i in xrange(N):
		if args.verbose:
			print('Calculating progress: '+str(int(float(i+1)*100/N))+'%',end="\r")
			sys.stdout.flush()
		accessibility[i] = calculateNodeDiversity(g,args.walk,i,args.steps,directed)
	if args.verbose:
		print("")


if args.outputfile:
	if args.verbose:
		print("Writing network to file \""+args.outputfile+"\"... ",end="")
	for i in xrange(N):
		g.vs[i]["accessibility"] = accessibility[i]	
	Graph.write(g,args.outputfile)
	if args.verbose:
		print("Done")
else:
	for i in xrange(N):
		print(str(i)+"\t"+str(accessibility[i]))
		
if args.verbose:
	print("Calculation done!")

# ----------------------------------------------------------------

