#!/bin/sh

# All values
in="gwas.rob_genes.simplify.txt"

cut -f 4 $in > bf.txt
paste bf.txt $in | sort -r -g | cut -f 2- > gwas.rob_genes.sort_bf.txt

cut -f 7 $in > pvalue.txt
paste pvalue.txt $in | sort -g | cut -f 2- > gwas.rob_genes.sort_pval.txt

# Same, using theta filtered values
in="gwas.rob_genes.simplify.theta.txt"

cut -f 4 $in > bf.txt
paste bf.txt $in | sort -r -g | cut -f 2- > gwas.rob_genes.sort_bf_theta.txt

cut -f 7 $in > pvalue.txt
paste pvalue.txt $in | sort -g | cut -f 2- > gwas.rob_genes.sort_pval_theta.txt
