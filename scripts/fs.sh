#!/bin/sh

for f in fs.txt.gz 1a17.txt.gz 1ubi.txt.gz
do
	out=`basename $f .txt.gz`_R.txt
	echo "IN: $f\t\tOUT: $out"

	echo "pdb.id\tchain.id\tdist\taa1\taa1.pos\taa2\taa2.pos\tmi\th\tvi\thxy\thyx\thx\thy\tconsx\tconsy\tcorr.naive\tcorr\tcorr.fodor\tllmsa" > $out
	gunzip -c $f | cut -f 1-7,19-30,33 >> $out
done

