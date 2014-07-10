#!/bin/sh

LIB=$HOME/workspace/Epistasis/lib/

java -Xmx10G -Djava.library.path=$LIB -jar $HOME/snpEff/Epistasis.jar $*
