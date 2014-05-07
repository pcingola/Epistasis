#!/bin/sh

numBases="3"
distance="3"

# # Perform PDB distance analysis
#java -jar Epistasis.jar pdbdist 8 db/pdb
java -jar Epistasis.jar pdbdist $distance db/pdb db/pfam/ensemblId_to_proteinId.txt  > pdbdist.$distance.txt
mv pdb_distance_by_AA_pos.txt pdb_distance_by_AA_pos.$distance.txt

# # Mutual information analysis
# java -jar Epistasis.jar \
# 	mi \
# 	100 \
# 	$numBases \
# 	/home/pcingola/snpEff/db/multiz100way/head.fa \
# 	2> epistasis.mi.$numBases.err \
# 	> epistasis.mi.$numBases.out

