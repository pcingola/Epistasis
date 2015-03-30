#!/usr/bin/perl

# Debug?
$debug = 0;

# Parse input
while( $l =<STDIN> ) {
	chomp $l;

	# Ignore comments
	if( $l !~ /^#/ ) {
		($pos1, $gene1, $ann1, $pos2, $gene2, $ann2, $t, $log10BFlogReg, $pvalue, $logLikCoEvo, $log10BFtotal) = split /\t/, $l;
    	print "($pos1, $gene1, $ann1, $pos2, $gene2, $ann2, $log10BFlogReg, $pvalue, $logLikCoEvo, $log10BFtotal)\n" if $debug;
	
		$score = $log10BFtotal;
		$score = -log($pvalue);

		if( $genes{$gene1} eq '' ) {
			print "$gene1\t$score\n";
			$genes{$gene1} = $score;
    	}
	
		if( $genes{$gene2} eq '' ) {
			print "$gene2\t$score\n";
			$genes{$gene2} = $score;
		}
	}
}
