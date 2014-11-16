#!/bin/sh

for i in */geneGeneTopValue_threshold*.txt
do
	base=`basename $i .txt`
	dir=`dirname $i`
	out="$dir/$base.tsv"

	echo $out
	cat $i | ./geneGeneTopValueThreshold2Txt.py > $out
done
