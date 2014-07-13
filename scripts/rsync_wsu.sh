#!/bin/sh

rsync -avz -e ssh eq8302@ehs.grid.wayne.edu:snpEff/epistasis/*.txt $HOME/snpEff/epistasis/
