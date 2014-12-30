#!/bin/sh

cat gwas.*.txt | ./gwasCoeff_altNull.pl \
	| tr ":][," "    " \
	| tr " " "\t" \
	| tr -s "\t" \
	| grep "Alt" \
	| cut -f 1,3- \
	> alt.txt

cat gwas.*.txt | ./gwasCoeff_altNull.pl \
	| tr ":][," "    " \
	| tr " " "\t" \
	| tr -s "\t" \
	| grep "Null" \
	| cut -f 1,3- \
	> null.txt
