#!/usr/bin/perl

use strict;

my( %genes );

#-------------------------------------------------------------------------------
# Is any of the genes in 'genesStr' resent in 'genes' hash?
#-------------------------------------------------------------------------------
sub hasGenes($) {
	my($genesStr) = @_;

	my($g);
	my(@t) = split /,/, $genesStr;
	foreach $g ( @t ) {
		if( $genes{$g} ne '' )	{ return 1; }
	}

	return 0;
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

# Read genes
my $genesFile = $ARGV[0];

# Read genes
die "Usage: cat gwas.20.BF.sort.ann.txt | ./filterGeneSet.pl genes_file.txt" if $genesFile eq '';

open GENES, $genesFile;
my $l;
while( $l = <GENES> ) {
	chomp $l;
	$genes{$l} = 1;
}
close GENES;

# Process STDIN
while( $l = <STDIN> ) {
	chomp $l;
	my($bf, $pval, $lltot, $lllr, $llmsa, $id1, $id2, $genes1, $ann1, $genes2, $ann2, $genesShared, $anns1, $anns2 ) = split /\t/, $l;

	# Filter lines
	if( hasGenes($genes1) && hasGenes($genes2) ) { print "$l\n"; }
	
}
