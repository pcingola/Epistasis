#!/usr/bin/perl

$absTh = 3.0;
$debug = 0;

while( $l =<STDIN> ) {
	chomp $l;
	@t = split /\t/, $l;

	# Parse log(Bf)
	($logBf, $logBfLogReg) = ( $t[3], $t[5] );

	# Show only if BF is increased by LL(MSA)
	$deltaLogBf = $logBf - $logBfLogReg;

	if( $deltaLogBf > 0 ) {
		# Parse logistic regression coefficients
		$thetaAlt = $t[19];
		$thetaAlt =~ tr/\[\]//d;
		@ths = split /,/, $thetaAlt;
		$sumAbs = abs($ths[0]) + abs($ths[1]);

		if( ( $absTh * $sumAbs ) < abs($ths[2]) ) {
			if( $debug )	{ print "$deltaLogBf\t$ths[0]\t$ths[1]\t$ths[2]\t$l\n"; } 
			else			{ print "$l\n"; }
		}
	}
}
