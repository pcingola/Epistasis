#!/bin/sh

# # Get all non-zero BF entries
# cat gwas.30.* | grep -v "log10(BF): 0.0" > gwas.logBf_non_zero.txt
# 
# # Extract theta parameters
# cat gwas.logBf_non_zero.txt | cut -f 18- > gwas.logBf_non_zero.theta.txt

# Split into ALT/NULL
cat gwas.logBf_non_zero.theta.txt | cut -f 1 | tr ":,][" "    " | tr -s " " | tr " " "\t" | cut -f 2- | sed -E "s/	$//g" > gwas.logBf_non_zero.theta.alt.txt
cat gwas.logBf_non_zero.theta.txt | cut -f 2 | tr ":,][" "    " | tr -s " " | tr " " "\t" | cut -f 2- | sed -E "s/	$//g" > gwas.logBf_non_zero.theta.null.txt
