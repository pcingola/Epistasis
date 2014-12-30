#!/usr/bin/perl

while( $l = <STDIN> ) {
	chomp $l;

	if( $l =~ /^\s+Alt.*:.*(\[.*\])/ ) {
		$alt = $1;
	} elsif( $l =~ /^\s+Null.*:.*(\[.*\])/ ) {
		$null = $1;
		print "$main\tLogReg_Alt: $alt\tLogReg_Null: $null\n";
	} else {
		$main = $l;
	}



	$prev = $l;
}
