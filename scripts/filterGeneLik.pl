#!/usr/bin/perl

#-------------------------------------------------------------------------------
# Filter out lines having low log-likelihood
#
#
# 															Pablo Cingolani 2014
#-------------------------------------------------------------------------------

$minll = 1.0;

while( $l = <STDIN> ) {
	($id1, $id2, $ll ) = split /\t/, $l;

	if( $ll > $minll )	{ print $l; }
}
