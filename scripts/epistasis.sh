#!/bin/sh

LIB=$HOME/snpEff/epistasis/lib/

java -Xmx10G -Djava.library.path=$LIB -jar $HOME/snpEff/Epistasis.jar $*
