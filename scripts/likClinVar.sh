#!/bin/sh

echo -e "ID\tll.max\tll.same"
for t in *txt
do
	cat $t | ./likClinVar.py > /dev/null
done
