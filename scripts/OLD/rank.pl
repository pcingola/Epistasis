#!/usr/bin/perl

while( $l = <STDIN> ) {
	chomp $l;
	($bf, $pval, $ll, $lllr, $llmsa, $id1, $id2, $genes1, $ann1, $genes2, $ann2) = split /\t/, $l;

	$score = $bf;
	$score = $llmsa;

	# Get a list of all genes in genes1 and genes2
	$genes12 = "$genes1\t$genes2";
	$genes12 =~ tr/,/\t/;
	@genes = split /\t/, $genes12;

	# Show score (if not done before)
	foreach $gene ( @genes ) {
		if( $scores{$gene} eq '' ) { 
			print "$gene\t$score\n"; 
			$scores{$gene} = $score;
		} 
	}

}
