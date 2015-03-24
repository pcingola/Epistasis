#!/usr/bin/perl

while( $l = <STDIN> ) {
	chomp $l;
	@t = split /\t/, $l;

	$talt = $t[19];
	$talt =~ tr/\[\]//d;
	@theta = split /,/, $talt;
	print "$theta[0] | $theta[1] | $theta[2]\t$talt\t$l\n" if $debug;

	$sumAbs = abs($theta[0]) + abs($theta[1]);
	$th3abs = abs($theta[2]);

	if( $th3abs > 3 * $sumAbs ) {
		print "[$theta[0] | $theta[1] | $theta[2] ]\t$l\n" if $debug;
		print "$l\n";
	}
}
