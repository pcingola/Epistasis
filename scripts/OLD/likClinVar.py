#!/usr/bin/env python

import sys

# Split an MSA id
def msaIdParse(id):
	id = id.replace(":","\t").replace("-","\t").replace("[","\t").replace("]","\t")
	(trid, start, end, aaidx, t) = id.split("\t")
	return (trid, start, end, aaidx)

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

linesByvcfId = dict()

# Read stdin
count = 1
llrsame = 0.0
llrmax = 0.0
msaId = ''
for l in sys.stdin:
	# Parse input
	(vcfId, msa1, msa2, llr, llnull, llalt) = l.rstrip().split('\t')
	(trid1, start1, end1, aaidx1) = msaIdParse(msa1)
	(trid2, start2, end2, aaidx2) = msaIdParse(msa2)

	ll = float(llr)

	if  msa1 == msa2: 
		llrsame = ll
		llr = "-1"
		msaId = msa1
	elif ll < 0: 
		llr = "0"
	else:
		llrmax = max(llrmax, ll)
		

	# Add to map
	print "\t".join( (trid1, start1, end1, aaidx1, trid2, start2, end2, aaidx2, llr, llnull, llalt) )

	count += 1


print >> sys.stderr, msaId + '\t' + str(llrmax) + '\t' + str(llrsame)
