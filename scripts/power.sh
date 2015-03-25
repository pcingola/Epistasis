#!/bin/sh

echo "lines n af1 af2 beta3 pmin p80 pmax" | tr " " "\t"

(
for f in n_*.txt
do
	cat $f | grep p-value | cut -f 7 | cut -f 2 -d : | tr -d " " | sort -g > tmp
	echo `cat tmp | wc -l` $f `head -n 1 tmp` `head -n 800 tmp | tail -n 1` `tail -n 1 tmp` 
done
) \
	| tr " " "\t" \
	| sed "s/.af1_/	/" \
	| sed "s/.af2_/	/" \
	| sed "s/.beta3_/	/" \
	| sed "s/.txt//" \
	| sed "s/n_//" \
