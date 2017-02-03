#!/usr/bin/python

from __future__ import print_function
import sys,argparse
from igraph import *
from Bio import SeqIO
from Bio import pairwise2

class Gene:
	
	def __init__(self, sequence, genDBname,annotation="",stringDBname=""):
		self.sequence = sequence
		self.genDBname = genDBname
		self.stringDBname = stringDBname
		self.annotation = annotation
		self.match = False
		
	def matchSequences(self,seq,name,concordance=0.9):
		alignScore = pairwise2.align.globalxx(self.sequence, seq, score_only = True)
		if len(self.sequence)>=len(seq):
			alignScore/=len(self.sequence)
		else:
			alignScore/=len(seq)
		if alignScore>=concordance:
			self.stringDBname = name
			self.match = True
			return True
		else:
			return False
	
			
class Genes:
	
	
	def __init__(self):
		self.geneList = []
		self.totalMatches = 0
		self.namesDict = {}
		self.annotationDict = {}
		
	def createGenesFromGB(self,gbfile):
		gb = SeqIO.parse(gbfile, "genbank")
		for seq_record in gb:
			for feat in seq_record.features:
				if 'old_locus_tag' in feat.qualifiers.keys() and 'translation' in feat.qualifiers.keys():
					gene = Gene(feat.qualifiers['translation'][0],feat.qualifiers['old_locus_tag'][0],annotation=feat.qualifiers['product'][0])
					self.geneList.append(gene)
		return self.geneList


	def getStringDBNames(self,stringfile,concordance=0.7):
		count = 0
		for seq_record in SeqIO.parse(stringfile, "fasta"):
			count+=1
			print('Matching progress: '+str(int(float(count+1)*100/len(self.geneList)))+'%',end='\r')
			sys.stdout.flush()
			for gene in self.geneList:
				if gene.match:
					continue
				if gene.matchSequences(seq_record.seq,seq_record.id,concordance):
					self.totalMatches += 1
					break
		print("")
		return self.totalMatches
	
	def generateDictionary(self):
		for gene in self.geneList:
			if gene.match:
				self.namesDict[gene.stringDBname] = gene.genDBname
		return self.namesDict
	
	def generateAnnotationDictionary(self):
		for gene in self.geneList:
			if gene.match:
				self.annotationDict[gene.genDBname] = gene.annotation
		return self.annotationDict

class Network:

	def __init__(self,genes):
		self.g = Graph()
		self.genes = genes
		
	def readStringNetworkFile(self,stringfile,measureCol="combined_score"):
		self.genes.generateDictionary()
		self.genes.generateAnnotationDictionary()
		f = open(stringfile)
		firstLine = True
		measureIndex = 0
		neighborhoodIndex = 0
		tuples = []
		for line in f.readlines():
			line = line.replace("\n","")
			line = line.split(" ")
			if firstLine:
				firstLine=False
				measureIndex = line.index(measureCol)
				neighborhoodIndex = line.index("neighborhood")
				continue
			measure = int(line[measureIndex])
			if measureCol == "combined_score":
				measure = int(line[measureIndex])-int(line[neighborhoodIndex])
			if line[0] in self.genes.namesDict.keys() and line[1] in self.genes.namesDict.keys() and measure:
				tuples.append([self.genes.namesDict[line[0]],self.genes.namesDict[line[1]],measure])
		self.g = Graph.TupleList(tuples,edge_attrs=[measureCol])
		for i in xrange(len(self.g.vs)):
			print(self.genes.annotationDict[self.g.vs[i]['name']])
			self.g.vs[i]['bdName'] = self.genes.annotationDict[self.g.vs[i]['name']]
		return self.g
		
	def writeNetworkToFile(self,filename):
		Graph.write(self.g,filename)
			
			

#Reading and treating arguments ----------------------------------

import argparse
parser = argparse.ArgumentParser(description="Calculate the accessibility value of each node of a network given as input.")
parser.add_argument("inputgb", help="Genbank file with genes and protein sequences")
parser.add_argument("inputstringlinks", help="String-db file with links")
parser.add_argument("inputstringseq", help="String-db file with protein sequences")
parser.add_argument("outputfile", help="File to output the generated network.")
parser.add_argument("-m","--measure", help="Measure from String-DB to consider as edge attribute. Default=combined_score",type=str,default="combined_score")
parser.add_argument("-v","--verbose", help="Make the software more verbose.",action="store_true")
args = parser.parse_args()



# ----------------------------------------------------------------

genes = Genes()
genes.createGenesFromGB(args.inputgb)
print(len(genes.geneList))

genes.getStringDBNames(args.inputstringseq,0.7)
print(genes.totalMatches)

net = Network(genes)
net.readStringNetworkFile(args.inputstringlinks,args.measure)
net.writeNetworkToFile(args.outputfile)

# ----------------------------------------------------------------

