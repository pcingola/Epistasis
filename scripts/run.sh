#!/bin/sh

echo Rmoving old files
rm nohup.out
rm -rvf epistasis.bds.* run.bds.* 

echo Process is run as nohup
nohup `dirname $0`/run.bds &

echo
echo
echo
