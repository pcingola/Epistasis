#!/usr/bin/env python

"""
  Select GeneSets from MSigDb using expression data from GTEx


  Pablo Cingolani 2013
"""

import sys

# Debug mode?
debug = False

# Value threshold
minMatchValue = float("-inf")
maxMatchValue = float("inf")
biogrid = dict()
reactome = dict()
geneId2Name = dict()

# Max number of missing IDs
maxMissingIds = 10

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
				if g1 in biogrid:
					biogrid[g1].add(g2)
				else:
					biogrid[g1] = set( [g2] )

				if debug: print >> sys.stderr, "BIOGRID: ", g1, '\t', g2, "\t", biogrid[g1]

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
		if debug : print >> sys.stderr, geneSetName, " => ", geneSet[ geneSetName ]
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
						if g1 in reactome:
							reactome[g1].add(g2)
						else:
							reactome[g1] = set( [g2] )

#------------------------------------------------------------------------------
# Process normalized GTEx file
#------------------------------------------------------------------------------
def readGtex(gtexFile, minMatchCount, minMatchValue, maxMatchValue, minAvgValue, maxAvgValue):
	print >> sys.stderr, "Reading GTEx file '{}'\n\tmin match: {}\n\tmin value: {}\n\tmax value: {}\n\tmin avg  : {}\n\tmax avg  : {}".format( gtexFile, minMatchCount, minMatchValue, maxMatchValue, minAvgValue, maxAvgValue )

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
			if missing > maxMissingIds: sys.exit(1)
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
	print >> sys.stderr, "Usage: " + sys.argv[0] + " msigDb.gmt gtex_normalized.txt gtexExperimentId_1,gtexExperimentId_2,...,gtexExperimentId_N minCount minMatchValue maxMatchValue minGeneSetSize"
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
minMatchCount = int( sys.argv[argNum] )

argNum += 1
minMatchValue = float( sys.argv[argNum] ) # Can be '-inf'

argNum += 1
maxMatchValue = float( sys.argv[argNum] ) # Can be 'inf'

argNum += 1
minAvgValue = float( sys.argv[argNum] ) # Can be '-inf'

argNum += 1
maxAvgValue = float( sys.argv[argNum] ) # Can be 'inf'

#---
# Load data
#---
loadGeneIds(geneId2NameFile)
loadReactomInt(reactomeIntFile)
loadBioGrid(bioGridFile)

# Parse IDs
ids = set( id for id in gtexExperimentIds.split(',') if id )	# Filter out empty IDs

# Read normalized GTEx file
gtexGenes = readGtex(gtexFile, minMatchCount, minMatchValue, maxMatchValue, minAvgValue, maxAvgValue)

# Select interactions for genes passing all filters
interactions = set()
for g1 in biogrid :
	if (g1 in reactome) and (g1 in gtexGenes):
		for g2 in biogrid[g1]:
			if (g2 in reactome) and (g2 in gtexGenes):
				if debug: print >> sys.stderr, "PASS:\t{}\t{}".format(g1, g2)

				if g1 == g2:
					print >> sys.stderr, "Skipping:\t{}\t{}".format(g1, g2)
				elif g1 < g2:
					interactions.add( "{}\t{}".format(g1, g2) )
				else:
					interactions.add( "{}\t{}".format(g2, g1) )

# Show results
for line in sorted( interactions ):
	print line
print >> sys.stderr, "Total number of pairs:", len( interactions )
