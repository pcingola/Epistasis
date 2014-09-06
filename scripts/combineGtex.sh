#!/bin/sh

DIR=$HOME/snpEff/epistasis
DIR_SCRIPTS=$HOME/snpEff/epistasis/scripts
DIR_GTEX=$HOME/snpEff/db/GRCh37/GTEx
DIR_REACTOME=$HOME/snpEff/db/reactome
DIR_BIOGRID=$HOME/snpEff/db/biogrid

# # Create geneID-geneName map file
# zcat $HOME/snpEff/data/GRCh37.75/genes.gtf.gz \
# 	| cut -f 3- \
# 	| grep "^gene" \
# 	| cut -f 7 \
# 	| tr \" "\t" \
# 	| cut -f 2,4 \
# 	| sort \
# 	| uniq \
# 	> geneId_geneName.txt

# # Select interations from Reactome (only direct_complex interactions)
# cat $DIR_REACTOME/interactions/homo_sapiens.interactions.txt \
# 	| grep -P "\tdirect_complex" \
# 	| cut -f 2,5 \
# 	| sed "s/ENSEMBL://g" \
# 	| sort \
# 	| uniq \
# 	> $DIR_REACTOME/reactome.homo_sapiens.interactions.txt

# # Get ID -> tissue mapping
# cut -f 1,7 $DIR_GTEX/gtex_ids.txt \
# 	| tr "-" "." \
# 	| sed "s/ . / - /" \
# 	> $DIR_GTEX/gtex_tissue.txt 

# # Select GTEx IDs that are related to pancreas
# cat $DIR_GTEX/gtex_tissue.txt \
# 	| grep Pancreas \
# 	| cut -f 1 \
# 	| tr "\n" "," \
# 	> $DIR_GTEX/pancreas_ids.txt

# Combine GTEX + Reactome + BioGrid
$DIR_SCRIPTS/combineGtex.py \
	$DIR/geneId_geneName.txt \
	$DIR_REACTOME/reactome.homo_sapiens.interactions.txt \
	$DIR_BIOGRID/biogrid.human.uniq.txt \
	$DIR_GTEX/gtex_norm.txt \
	`cat $DIR_GTEX/pancreas_ids.txt` \
	5 \
	0 \
	inf \
	0.3 \
	inf \

	
