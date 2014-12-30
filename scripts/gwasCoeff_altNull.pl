#!/usr/bin/perl


$debug = 0;
$showNext = 0;

while( $l = <STDIN> ) {
	chomp $l;

	if( $l =~ /log\(BF_LogReg\): (\S+)/ ) {
		$lbf = $1;
		$showNext = 0;
		if( $lbf > 0)	{ 
			$showNext = 1;
			print "$l\n\t$lbf\n" if $debug; 
		}
	} elsif( $showNext ) {
		print "$lbf\t$l\n";
	}

	
}
