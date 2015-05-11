#!/bin/sh -e

# Create thesis
( pdflatex thesis_mcgill ; bibtex thesis_mcgill) 2>&1 | tee make.out

echo
echo "================================================================================"
echo " Citations"
echo "================================================================================"
grep Citation make.out || true

echo
echo "================================================================================"
echo " References"
echo "================================================================================"
grep Reference make.out || true

open *.pdf

