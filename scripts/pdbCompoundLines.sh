#!/bin/sh

for pdb in `find . -iname "*.pdb"`
do
	pdbId=`basename $pdb .pdb`
	cat $pdb | /home/pcingola/workspace/Epistasis/scripts/pdbCompoundLines.py $pdbId
done
