#!/bin/sh

echo Rmoving old files
rm nohup.out
rm -rvf epistasis.bds.*

echo Process is run as nohup
nohup ./scripts/run.bds &

#echo Tailing nohup.out, you can Ctrl-C
#sleep 1
#tail -f nohup.out 

