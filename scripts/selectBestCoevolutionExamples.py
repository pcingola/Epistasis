#!/usr/bin/env python

#-------------------------------------------------------------------------------
#
# Select co-evolution examples from a 'likelihood.pdb_compound.*.txt' file
#
#
#															Pablo Cingolani 2015
#-------------------------------------------------------------------------------

import sys

# A set consisting of a GAP character
setGap = set('-')

# Do all bases in the sequence have the same base (gaps are ignored)
def sameBase(seq):
	charSet = set(seq) - setGap 
	return len(charSet) <= 1

# Count pairs of charates (one form each sequence)
# Ignore pair if either AA is a gap
def pairs(seq1, seq2):
	count = dict()

	for pair in zip(seq1, seq2):
		# Neither are gaps?
		if pair[0] != '-' and pair[1] != '-':
			# Join them to get 'AA-pair'
			aa = "".join( sorted(pair) )

			# Count entries
			if aa in count: count[aa] += 1
			else:			count[aa] = 1

	# Sort entries by count value
	swapCounts = [ (count[k], k) for k in count.keys() ]
	countSort = sorted(swapCounts, reverse=True)

	# Show
	out = ""
	for d in countSort:
		out += '\t\t' + str(d) + '\n'

	# Do we satisfy filter condition?
	if len(countSort) > 1:
		first = countSort[0]
		second = countSort[1]

		# Does the second entry have than one item?
		countSecond = second[0]
		if countSecond > 1:
			# Do top-2 entries refer to different AAs?
			if len(set(first[1] + second[1])) == 4:
				return (True, out)

	return (False, out)

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------
for line in sys.stdin:
	(dist, aa1, aa2, llr, llint, llnull, seq1, seq2, etc) = line.split('\t', 8)

	llr = float(llr)
	if llr > 1 and not sameBase(seq1) and not sameBase(seq2):
		(ok, out) = pairs(seq1, seq2)
		if ok :
			print line.rstrip()
			print llr, '\n\tseq1: ', seq1, '\n\tseq2: ', seq2
			print '\tAA pairs:\n', out, "\tOK:", ok, '\n'

