#!/bin/sh

echo "dist\tmi\th\tvi\thxy\thyx\thx\thy\tconsx\tconsy\tcorr.naive\tcorr\tcorr.fodor\tllmsa"
gunzip -c fs.txt.gz | cut -f 3,19-30,33
