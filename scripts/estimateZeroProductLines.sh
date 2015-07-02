#!/bin/sh

head -n 100000 gwas.alleleMat.txt | ./estimateZeroProductLines.pl
