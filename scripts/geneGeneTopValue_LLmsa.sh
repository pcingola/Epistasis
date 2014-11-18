#!/bin/sh

for txt in `find . -iname "*.txt"`
do
	max=`cat $txt | ./geneGeneTopValue_LLmsa.pl`
	echo $txt $max | tr "/ " "\t\t" | sed "s/.txt//" | cut -f 2- | tee -a geneGeneTopValue_LLmsa_2.txt
done
