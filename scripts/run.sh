#!/bin/sh

echo Rmoving old files
rm nohup.out
rm -rvf epistasis.bds.* run.bds.* 

echo Process is run as nohup
nohup ./scripts/run.bds &

echo
echo
echo
