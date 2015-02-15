#!/usr/bin/env python

import sys

print '\t'.join( ('gene1', 'gene2', 'threshold_ll_null', 'll_max', 'count', 'countPass', 'id1', 'id2', 'll', 'llalt', 'llnull', 'seq1', 'seq2') )
for l in sys.stdin:
	f = l.replace('/', '\t').replace(' ','\t').replace('.txt','').rstrip().split('\t')
	if len(f) <= 17: continue

	(gene1, gene2, threshold_ll_null, ll_max, count, countPass, id1, id2, ll, llalt, llnull, seq1, seq2) = (f[1], f[2], f[4], f[6], f[8], f[10], f[11], f[12], f[13], f[14], f[15], f[16], f[17])
	print '\t'.join( (gene1, gene2, threshold_ll_null, ll_max, count, countPass, id1, id2, ll, llalt, llnull, seq1, seq2) )
