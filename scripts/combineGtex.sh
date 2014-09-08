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

# # Select interactios from BioGrid
# # Note: We filter out RPS/RPL robosomal complexes interactions
# cat $DIR_BIOGRID/BIOGRID-ALL-3.2.115.tab2.txt \
# 	| cut -f 8,9,16,17 \
# 	| grep -P "\t9606\t9606" \
# 	| cut -f 1,2 \
# 	| sort \
# 	| uniq \
# 	| grep -v "^RPS" \
# 	| grep -v "^RPL" \
# 	| grep -vP "\tRPS" \
# 	| grep -vP "\tRPL" \
# 	> $DIR_BIOGRID/biogrid.human.uniq.txt

# # Get ID -> tissue mapping
# cut -f 1,7 $DIR_GTEX/gtex_ids.txt \
# 	| tr "-" "." \
# 	| sed "s/ . / - /" \
# 	> $DIR_GTEX/gtex_tissue.txt 

# # Gene names in MSAs
# cut -f 3 $DIR/idMap_ensemblId_refseq_pdbId.best.txt | sort | uniq > $DIR/genes.msas.txt

# Calculate combined interaction pairs
rm -rvf $DIR/interactions.*.txt
for TISSUE in "Adipose - Subcutaneous" "Adipose - Visceral" "Liver" "Muscle - Skeletal" "Pancreas"
do

	TIS=`echo $TISSUE | tr -d "\n-" | tr " " "_" | tr -s "_"`
	echo Tissue: $TISSUE $TIS

	# Select GTEx IDs that are related to pancreas
	cat $DIR_GTEX/gtex_tissue.txt \
		| grep "$TISSUE" \
		| cut -f 1 \
		| tr "\n" "," \
		> $DIR_GTEX/$TIS.ids.txt

	# Combine GTEX + Reactome + BioGrid
	$DIR_SCRIPTS/combineGtex.py \
		$DIR/geneId_geneName.txt \
		$DIR_REACTOME/reactome.homo_sapiens.interactions.txt \
		$DIR_BIOGRID/biogrid.human.uniq.txt \
		$DIR/genes.msas.txt \
		$DIR_GTEX/gtex_norm.txt \
		`cat $DIR_GTEX/$TIS.ids.txt` \
		0.25 \
		0 \
		inf \
		0.3 \
		inf \
		| tee $DIR/interactions.$TIS.txt

done

cat $DIR/interactions.*.txt | sort | uniq > $DIR/interactions.txt

