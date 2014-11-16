#!/bin/sh

for th in 1E-100 1E-80 1E-60 1E-50 1E-40 1E-30 1E-20 1E-10 
do
	for txt in `find . -iname "*.txt"`
	do
		max=`cat $txt | ./geneGeneTopValue_threshold.pl $th`
		echo $txt $max | tee -a geneGeneTopValue_threshold.$th.txt
	done
done
