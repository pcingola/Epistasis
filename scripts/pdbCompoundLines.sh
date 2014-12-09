#!/bin/sh

dir=$1

for pdb in `find $dir -iname "*.pdb"`
do
	pdbId=`basename $pdb .pdb`
	cat $pdb | $HOME/snpEff/epistasis/scripts/pdbCompoundLines.py $pdbId
done
