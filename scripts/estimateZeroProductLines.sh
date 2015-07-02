#!/bin/sh

head -n 10000 gwas.alleleMat.txt | ./estimateZeroProductLines.pl
