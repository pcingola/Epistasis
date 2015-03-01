#!/bin/sh

# Delete old 'tmp' files
#rm -vf *.bbl *.aux *.blg *.log *.lof *.lot *.pdf *.synctex.gz

( pdflatex thesis_mcgill.tex ; bibtex thesis_mcgill) 2>&1 | tee make.out

open thesis_mcgill.pdf

