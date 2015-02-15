#!/bin/sh

rm -vf top.txt

for t in by_variant/*.txt
do
	var=`head -n 1 $t | cut -f 1 | tr ":_/" "\t\t\t"`
	top=`cut -f 4 $t | sort -rg | head -n 1`
	wcl=`wc -l $t`

	echo $var $top $wcl | tee -a top.txt
done

cat top.txt | tr " " "\t" | cut -f 1-6 | sort -k 1n,1n -k 2n,2n  > top.sorted.txt
