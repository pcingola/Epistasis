#!/bin/sh

#-------------------------------------------------------------------------------
# Note: Most of these files are calculated in scripts/combineGtex.sh
#-------------------------------------------------------------------------------

./scripts/selectNotInteracting.py \
	geneId_geneName.txt \
	reactome.homo_sapiens.interactions.txt \
	biogrid.human.all.txt \
	c2.cp.v4.0.symbols.gmt \
	genes.msas.txt \
	| head -n 1000000 \
	| tee selectNotInteracting.txt

cat selectNotInteracting.txt | grep "^OK" | cut -f 2,3 | head -n 10000 > selectNotInteracting.head.txt 
