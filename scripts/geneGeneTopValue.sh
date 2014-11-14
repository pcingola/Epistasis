#!/bin/sh

for txt in `find likelihood_genes -iname "*.txt"`
do
	max=`cut -f 3 $txt | ./topValue.pl`
	echo $txt $max | tee -a topValues.txt
done
