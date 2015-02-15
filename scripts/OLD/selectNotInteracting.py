#!/usr/bin/env python

"""
  Select pairs of genes that are assumed NOT to interact:
	- Not in PPI database (BioGrid)
	- Not in the same pathway (GeneSet: C2 from MSigDb)
	- Not in interacting in Reactme 

  Pablo Cingolani 2014
"""

import sys
import random

# Debug mode?
debug = False

# Value threshold
interactions = dict()		# Keep all interactions in this dictionary
geneId2Name = dict()	# Gene ID conversion map

# Max number of missing IDs
maxMissingIdsFactor = 0.5

#------------------------------------------------------------------------------
# Add interaciotn to hash
#------------------------------------------------------------------------------
def addInteraction(g1, g2, name):
	key = genes2key(g1, g2)

	if key not in interactions:
		interactions[key] = set()

	interactions[key].add(name)
		

#------------------------------------------------------------------------------
# Create a string 'key' from a gene pair
#------------------------------------------------------------------------------
def genes2key(g1, g2):
	if( g1 > g2 ): return g2 + "\t" + g1
	return g1 + "\t" + g2
	
#------------------------------------------------------------------------------
# Load gene ID <-> name mapping
#------------------------------------------------------------------------------
def loadBioGrid(biogridFile):
	print >> sys.stderr, "Loading BioGrid interactions from '{}'".format( biogridFile )

	for line in open(biogridFile) :
		fields = line.rstrip().split("\t")
		if len(fields) == 2:
			(g1, g2) = fields
			if g1 and g2: addInteraction(g1, g2, "BioGrid")

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
# Load gene set (GMT format from MSigDb)
#------------------------------------------------------------------------------
def loadGeneSets(geneSetsFile):
	print >> sys.stderr, "Loading GeneSets from file '{}'".format( geneSetsFile )

	for line in open(geneSetsFile) :
		fields = line.rstrip().split("\t")
		if len(fields) >= 3:
			geneSetName = fields[0]
			print "\t", geneSetName
			for n1 in range(2, len(fields)):
				for n2 in range(n1+1, len(fields)):
					g1 = fields[n1]
					g2 = fields[n2]
					addInteraction(g1, g2, geneSetName )

#------------------------------------------------------------------------------
# Load MSigDb file
#------------------------------------------------------------------------------
def loadMsigDb(msigFile):
	geneSet = {}
	for line in open(msigFile) :
		fields = line.rstrip().split("\t")
		geneSetName = fields[0]
		geneSet[ geneSetName ] = fields[2:]
		if debug : print >> sys.stderr, "MSIGDB:\t", geneSetName, " => ", geneSet[ geneSetName ]
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

						addInteraction(g1, g2, "Reactome_interaction")

#------------------------------------------------------------------------------
# Process normalized GTEx file
#------------------------------------------------------------------------------
def readGtex(gtexFile, minMatchPercent, minMatchValue, maxMatchValue, minAvgValue, maxAvgValue):
	minMatchCount = minMatchPercent * len(ids)
	print >> sys.stderr, "Reading GTEx file '{}'\n\tmin match: {}% ( {} )\n\tmin value: {}\n\tmax value: {}\n\tmin avg  : {}\n\tmax avg  : {}".format( gtexFile, 100 * minMatchPercent, minMatchCount, minMatchValue, maxMatchValue, minAvgValue, maxAvgValue )

	columnIdx = []
	header = []
	countOk = 0
	countNo = 0
	gtexGenes = set()
	for line in open(gtexFile) :
		fields = line.rstrip().split("\t")

		if not columnIdx : 
			header = fields[:]
			# Read header and add all columnIdx numbers that we are looking for
			foundIds = set()
			for i in range(len(fields)):
				if fields[i] in ids:
					columnIdx.append(i)
					foundIds.add( fields[i] )

			# Sanity check
			missing = 0
			for id in ids:
				if id not in foundIds:
					print >> sys.stderr, "Missing GTEx ID '{}'.".format( id )
					missing += 1

			# Stop if too many are missing
			print >> sys.stderr, 'Missing samples: ', missing, ' / ', len(ids)
			if missing > (maxMissingIdsFactor * len(ids)): sys.exit(1)
			if debug: print >> sys.stderr, "OK, All required IDs found."

		else :
			geneId, geneName = fields[0], fields[1]

			# Collect values for requested IDs
			avg = 0
			count = 0
			missing = 0 
			vals = []
			for idx in columnIdx :
				val = fields[idx]
				if val != "NA":
					v = float(val)
					avg += v
					count += 1
					if (v >= minMatchValue) and (v <= maxMatchValue): vals.append( v )
				else:
					missing += 1

			# Show results
			if count > 0:	avg = avg / count
			if len(vals) >= minMatchCount and (minAvgValue <= avg) and (avg <= maxAvgValue):
				countOk += 1
				gtexGenes.add( geneName )
				if debug: print >> sys.stderr, "OK\t{}\t{}\tmatch: {}\tmissing: {}\ttotal: {}\tavg: {}\tvalues: {}".format(geneId, geneName, len(vals), missing, len(columnIdx), avg, vals)
			else:
				countNo += 1
				if debug: print >> sys.stderr, "NO\t{}\t{}\tmatch: {}\tmissing: {}\ttotal: {}\tavg: {}\tvalues: {}".format(geneId, geneName, len(vals), missing, len(columnIdx), avg, vals)

	print >> sys.stderr, "\tNumber of genes that passed  : {}\n\tNumber of genes filtered out : {}\n".format( countOk, countNo)
	return gtexGenes

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------

#---
# Command line parameters
#---
if len(sys.argv) < 4 :
	print >> sys.stderr, "Usage: " + sys.argv[0] + " geneId2NameFile.txt reactomeIntFile.txt bioGridFile.txt geneSet.MSigDb.gmt genesMsasFile.txt"
	sys.exit(1)

argNum=1
geneId2NameFile = sys.argv[argNum]
argNum += 1
reactomeIntFile = sys.argv[argNum]
argNum += 1
bioGridFile = sys.argv[argNum]
argNum += 1
geneSetsFile = sys.argv[argNum]
argNum += 1
genesMsasFile = sys.argv[argNum]

#---
# Load data
#---
loadGeneIds(geneId2NameFile)
loadReactomInt(reactomeIntFile)
loadBioGrid(bioGridFile)
loadGeneSets(geneSetsFile)
genesMsas = set( [ line.rstrip() for line in open(genesMsasFile) ] )
print >> sys.stderr, "Loaded genes (MSAs) from ", genesMsasFile, ". Genes loadded: ", len(genesMsas)

#---
# Select random pairs not in any set
#---

genes = list(genesMsas)
numGenes = len(genes)
while 1:
	n1 = random.randint(0, numGenes-1)
	g1 = genes[n1]
	n2 = random.randint(0, numGenes-1)
	g2 = genes[n2]

	key = genes2key(g1, g2)
	if key in interactions: print "NO\t", g1, "\t", g2, "\t", interactions[key]
	else: print "OK\t", g1, "\t", g2
	
