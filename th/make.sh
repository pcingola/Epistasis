#!/bin/sh -e

# Create thesis
( pdflatex thesis_mcgill ; bibtex thesis_mcgill) 2>&1 | tee make.out

echo
echo
echo
grep Citation make.out
grep Reference make.out

open *.pdf

