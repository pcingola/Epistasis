#!/usr/bin/perl

while( $l = <STDIN> ) {
	chomp $l;
	#($id1, $id2, $ll, $ln, $la, $s1, $s2) = split /\t/, $l;
	($ll, $ln, $la, $s1, $s2) = split /\t/, $l;
	print "$ll\t$ln\t$la\n\t$s1\n\t$s2\n\n";
}
