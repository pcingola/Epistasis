#!/bin/sh

distance="5.0"

java -Xmx4G -jar Epistasis.jar \
		checkPdbGenome \
		/home/pcingola/snpEff/snpEff.config  \
		hg19 \
		/home/pcingola/snpEff/db/pdb/pdb_hires_human/ \
		/home/pcingola/snpEff/db/multiz100way/hg19.100way.nh \
		/home/pcingola/snpEff/db/multiz100way/head.fa \
		/home/pcingola/snpEff/db/multiz100way/idMap_ensemblId_refseq_pdbId.txt \
		5.0 \
	> epistasis.out 2>&1

