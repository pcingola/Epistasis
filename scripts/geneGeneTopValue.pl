#!/usr/bin/perl

for( $count = 0 ; $l = <STDIN> ; $count++ ) {
	chomp $l;
	if( $count == 0 )	{ $max = $l; }
	elsif( $l > $max )	{ $max = $l; }
	
}

print "max: $max\tcount: $count\n";
