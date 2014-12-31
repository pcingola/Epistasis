#!/bin/sh

echo Adjusting Bayes Factors
# java -Xmx4G -jar $HOME/snpEff/Epistasis.jar adjustrawbayesfactor gwasCoeff_altNull.txt gwas.30.txt > gwas.30.adjusted.txt

echo Removing columns
# echo "#Log10(BF)\tLog10(BF_MSA)\tLog10(pThetaRatio)\tpValue\tID_i\tID_j\tAnn_i\tAnn_j" > gwas.main.txt
# cat gwas.30.adjusted.txt \
# 	| sed "s/ , / /g" \
# 	| tr " ][" "\t\t\t" \
# 	| tr -s "\t" \
# 	| cut -f 1,3,5,13,26- \
# 	>> gwas.main.txt

echo
echo Sorting by BF
#sort -rg gwas.main.txt | ./addGeneCol.pl > gwas.main.sort_bf.txt

echo
echo Sorting by p-value
#sort -k4g,4g gwas.main.txt | ./addGeneCol.pl > gwas.main.sort_pvalue.txt

echo
echo Top genes BF
cut -f 1,2 gwas.main.sort_bf.txt | head -n 100 | tr ",\t" "\n\n" | sort | uniq -c | sort -rn | head 

echo
echo Top genes p-value
cut -f 1,2 gwas.main.sort_pvalue.txt | head -n 100 | tr ",\t" "\n\n" | sort | uniq -c | sort -rn | head 

echo
echo Top variants BF
head -n 100 gwas.main.sort_bf.txt | cut -f 7,8 | tr "\t" "\n" | sort | uniq -c | sort -rn | head

echo
echo Top variants p-value
head -n 100 gwas.main.sort_pvalue.txt | cut -f 7,8 | tr "\t" "\n" | sort | uniq -c | sort -rn | head
