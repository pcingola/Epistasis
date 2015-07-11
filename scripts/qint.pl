#!/usr/bin/perl

$k = 10000;

while( $l = <STDIN> ) {
	chomp $l;
	@t = split /\t/, $l;
	foreach $n ( @t ) {
		if( $n =~ /_/ )	{ print "$n\t"; }
		else { 
			$n = int( $n * $k );
			if( $n == 0 )	{ print " \t"; } 
			else			{ print "$n\t"; }
		}
	}
	print "\n";
}
