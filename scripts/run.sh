#!/bin/sh

echo Rmoving old files
rm -rvf nohup.out epistasis.bds.* run.bds.* 

echo Process is run as nohup
nohup `dirname $0`/run.bds -model &

echo
echo
echo
