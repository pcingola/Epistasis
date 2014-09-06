#!/usr/bin/env python

"""
  Select GeneSets from MSigDb using expression data from GTEx


  Pablo Cingolani 2013
"""

import sys

# Debug mode?
debug = False

# Value threshold
minValue = float("-inf")
maxValue = float("inf")
interactions = dict()
geneId2Name = dict()

#------------------------------------------------------------------------------
# Load gene ID <-> name mapping
#------------------------------------------------------------------------------
def loadBioGrid(biogridFile):
	print >> sys.stderr, "Loading BioGrid interactions from '{}'".format( biogridFile )

	for line in open(biogridFile) :
		fields = line.rstrip().split("\t")
		if len(fields) == 2:
			(g1, g2) = fields

			# Convert gene IDs to gene names
			if g1 and g2:
				# Add to set
				if g1 in interactions:
					interactions[g1].add(g2)
				else:
					interactions[g1] = set(g2)

#------------------------------------------------------------------------------
# Load gene ID <-> name mapping
#------------------------------------------------------------------------------
def loadGeneIds(geneFile):
	print >> sys.stderr, "Loading gene ID <-> Name mapping '{}'".format( geneFile )

	for line in open(geneFile) :
		fields = line.rstrip().split("\t")
		if len(fields) == 2:
			(geneId, geneName) = fields
			geneId2Name[ geneId ] = geneName

#------------------------------------------------------------------------------
# Load MSigDb file
#------------------------------------------------------------------------------
def loadMsigDb(msigFile):
	geneSet = {}
	for line in open(msigFile) :
		fields = line.rstrip().split("\t")
		geneSetName = fields[0]
		geneSet[ geneSetName ] = fields[2:]
		if debug : print geneSetName, " => ", geneSet[ geneSetName ]
	return geneSet

#------------------------------------------------------------------------------
# Load Reactome interactions file
#------------------------------------------------------------------------------
def loadReactomInt(rintFile):
	print >> sys.stderr, "Loading interactions from '{}'".format( rintFile )

	for line in open(rintFile) :
		fields = line.rstrip().split("\t")
		if len(fields) == 2:
			(genes1, genes2) = fields

			# Multiple interactions are separated by '|'
			for g1 in genes1.split('|'):
				for g2 in genes2.split('|'):
					
					# Convert gene IDs to gene names
					if g1 and g2 and (g1 in geneId2Name) and (g2 in geneId2Name):
						g1 = geneId2Name[g1]
						g2 = geneId2Name[g2]

						# Add to set
						if g1 in interactions:
							interactions[g1].add(g2)
						else:
							interactions[g1] = set(g2)

#------------------------------------------------------------------------------
# Process normalized GTEx file
#------------------------------------------------------------------------------
def readGtex(gtexFile, minValCount, minValue, maxValue):
	print >> sys.stderr, "Reading GTEx file '{}'\n\tmin count: {}\n\tmin value: {}\n\tmax value: {}".format( gtexFile, minValCount, minValue, maxValue )

	columnIdx = []
	header = []
	countOk = 0
	countNo = 0
	gtexGenes = {}
	for line in open(gtexFile) :
		fields = line.rstrip().split("\t")

		if not columnIdx : 
			header = fields[:]
			# Read header and add all columnIdx numbers that we are looking for
			for i in range(len(fields)):
				if fields[i] in ids:
					columnIdx.append(i)
					ids[ fields[i] ] = 1

			# Sanity check
			ok = True
			for id in ids:
				if not ids[id]:
					print >> sys.stderr, "Missing GTEx ID '{}'.".format( id )
					ok = False

			if not ok: sys.exit(1)
			if debug: print >> sys.stderr, "OK, All required IDs found."

			# We require at least these number of values
			print >> sys.stderr, "Filter:\n\tMinimum number of values: {}\n\tMinimum value: {}\n\tMaximum value: {}\n".format(minValCount, minValue, maxValue)

		else :
			geneId, geneName = fields[0], fields[1]

			# Collect values for requested IDs
			vals = []
			for idx in columnIdx :
				val = fields[idx]
				if val != "NA":
					 v = float(val)
					 if (v >= minValue) and (v <= maxValue): vals.append( v )

			# Show results
			if len(vals) >= minValCount:
				avg = reduce(lambda x,y : x+y, vals) / len(vals)
				countOk += 1
				gtexGenes[geneName] = 1
				if debug: print "OK\t{}\t{}\t{} / {}\t{}\t{}".format(geneId, geneName, len(vals), len(columnIdx), avg, vals)
			else:
				countNo += 1
				gtexGenes[geneName] = 0
				if debug: print "NO\t{}\t{}\t{} / {}\t{}".format(geneId, geneName, len(vals), len(columnIdx), vals)

	print >> sys.stderr, "\tNumber of genes that passed  : {}\n\tNumber of genes filtered out : {}\n".format( countOk, countNo)
	return gtexGenes

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------

#---
# Command line parameters
#---
if len(sys.argv) < 4 :
	print >> sys.stderr, "Usage: " + sys.argv[0] + " msigDb.gmt gtex_normalized.txt gtexExperimentId_1,gtexExperimentId_2,...,gtexExperimentId_N minCount minValue maxValue minGeneSetSize"
	sys.exit(1)

argNum=1
geneId2NameFile = sys.argv[argNum]

argNum += 1
reactomeIntFile = sys.argv[argNum]

argNum += 1
bioGridFile = sys.argv[argNum]

argNum += 1
gtexFile = sys.argv[argNum]

argNum += 1
gtexExperimentIds = sys.argv[argNum]

argNum += 1
minValCount = int( sys.argv[argNum] )

argNum += 1
minValue = float( sys.argv[argNum] ) # Can be '-inf'

argNum += 1
maxValue = float( sys.argv[argNum] ) # Can be 'inf'

#---
# Load data
#---
loadGeneIds(geneId2NameFile)
loadReactomInt(reactomeIntFile)
loadBioGrid(reactomeIntFile)

gtexGenes = readGtex(gtexFile, minValCount, minValue, maxValue)

