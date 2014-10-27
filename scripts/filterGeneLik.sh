#!/bin/sh

rm likelihood_genes.txt

for g in */*.txt
do
	echo $g
	cat $g | ./filterGeneLik.pl >> likelihood_genes.txt
done
