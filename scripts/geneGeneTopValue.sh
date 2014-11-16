#!/bin/sh

for txt in `find . -iname "*.txt"`
do
	max=`cut -f 3 $txt | ./geneGeneTopValue.pl`
	echo $txt $max | tee -a topValues.txt
done
