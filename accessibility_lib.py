from igraph import *
import sys, getopt
import numpy as np

#Perform a self avoiding random walk in an igraph network net with S steps starting from the node i. Parameter directed = True if the network is directed
def selfAvoidWalk(net,S,i,directed=True):
	path = []
	node = net.vs[i]
	particle = i
	directed_mode = ""
	if directed:
		directed_mode = "OUT"
	else:
		directed_mode = "ALL"
	for n in xrange(S):
		neighbors = net.neighbors(net.vs[particle],mode=directed_mode)
		if i in neighbors:
			neighbors.remove(i)
		available = [e for e in neighbors if e not in path]
		if len(available):
			next = np.random.randint(len(available))
			path.append(available[next])
			particle = available[next]
		else:
			path.append(particle)
	return path
	
	
#Calculate the Diversity signature calculating diversity values from h=1..S from the node i in the net igraph network running M walks
def calculateNodeDiversitySignature(net,S,i,M=200,directed=True):
	N = len(net.vs)
	V = np.zeros(shape=(S,N))
	for n in xrange(M):
		path = selfAvoidWalk(net,S,i,directed)
		for h in xrange(S):
			V[h][path[h]] += 1
	diversity = np.zeros(shape=(S,1))	
	for h in xrange(S):
		nonZeroV = V[h][V[h]>0]
		probabilities = nonZeroV/M
		diversity[h] = -sum(probabilities*np.log(probabilities))
	return diversity

#Calculate the Diversity value from the node i in the net igraph network running M walks
def calculateNodeDiversity(net,h,i,M=200,directed=True):
	N = len(net.vs)
	V = np.zeros(shape=(N))
	for n in xrange(M):
		path = selfAvoidWalk(net,h,i,directed)
		V[path[-1]] += 1
	nonZeroV = V[V>0]
	probabilities = nonZeroV/M
	diversity = -sum(probabilities*np.log(probabilities))
	return diversity

#Calculate the Normalized Diversity value from the node i in the net igraph network running M walks
def calculateNormalizedNodeDiversity(net,h,i,M=200,directed=True):
	return calculateNodeDiversity(net,h,i,M,directed)/(1/np.log(len(net.vs)))



