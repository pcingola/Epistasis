#!/bin/sh

for pdb in `find . -iname "*.pdb"`
do
	pdbId=`basename $pdb .pdb`
	cat $pdb | $HOME/snpEff/epistasis/scripts/pdbCompoundLines.py $pdbId
done
