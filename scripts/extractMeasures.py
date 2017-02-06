#!/usr/bin/python

from __future__ import print_function
import sys,argparse
from igraph import *


#Reading and treating arguments ----------------------------------

import argparse
parser = argparse.ArgumentParser(description="Get measures from two networks and table the measures acording to the same nodes in both network.")
parser.add_argument("inputwcnet", help="Whole-cell network file")
parser.add_argument("inputstringnet", help="String-db network file")
parser.add_argument("outputfile", help="File to output the data.")
parser.add_argument("-v","--verbose", help="Make the software more verbose.",action="store_true")
args = parser.parse_args()



# ----------------------------------------------------------------

gWC = Graph.Load(args.inputwcnet)

gSDB = Graph.Load(args.inputstringnet)

gSDB.simplify()

gWC.vs["degree"] = gWC.degree()
gSDB.vs["degree"] = gSDB.degree()
gWC.vs["bt"] = gWC.betweenness()
gSDB.vs["bt"] = gSDB.betweenness()

prots = gWC.vs.select(annotation_in=["Protein Monomer"])


f = open(args.outputfile,"w")

f.write("NodeName\tWC_degre\tSDB_degree\tWC_btw\tSDB_btw\tWC_acs\tSDB_acs\n")

for prot in prots:
	stringNode = gSDB.vs.select(name_eq=prot['name'].replace("_MONOMER",""))
	if len(stringNode):
		stringNode = stringNode[0]
		f.write(stringNode['name']+"\t"+str(prot['degree'])+"\t"+str(stringNode['degree'])+"\t"+str(prot['bt'])+"\t"+str(stringNode['bt'])+"\t"+str(prot['accessibility'])+"\t"+str(stringNode['accessibility'])+"\n")


f.close()
