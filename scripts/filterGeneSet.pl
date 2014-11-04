#!/usr/bin/perl

#-------------------------------------------------------------------------------
#
# Filter a GWAS results using GeneSet
#
#
#
#															Pablo cingolani 2014
#-------------------------------------------------------------------------------

use strict;

my( %genes );

#-------------------------------------------------------------------------------
# Is any of the genes in 'genesStr' resent in 'genes' hash?
#-------------------------------------------------------------------------------
sub hasGenes($) {
	my($genesStr) = @_;

	# Split gene 'string' to get individual genes
	my($g, $s);
	my(@t) = split /,/, $genesStr;

	# Extract all gene sets from each gene
	my %sets = ();
	foreach $g ( @t ) {
		if( $genes{$g} ne '' )	{ 
			my(@gsets) = split /,/, $genes{$g};
			foreach $s ( @gsets ) { $sets{$s} = 1; } 
		}
	}

	# Join all 
	return join(',', sort keys %sets );
}

#-------------------------------------------------------------------------------
# Read GMT file
#-------------------------------------------------------------------------------
sub readGmt($) {
	my($gmt) = @_;

	open GENES, $gmt;

	my ($l, $i);
	while( $l = <GENES> ) {
		# Split by tab
		chomp $l;
		my(@genes) = split /\t/, $l;

		# First field is gene set name
		my($setName) = $genes[0];

		# Fields 3 and on are genes
		for( $i=2 ; $i <= $#genes ; $i++ ) {
			my($g) = $genes[$i];

			# Append gene set name
			if( $genes{$g} ne '' )	{ $genes{$g} .= ","; }
			$genes{ $g } .= $setName;
		}
	}

	close GENES;
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

# Read genes
my $genesFile = $ARGV[0];

# Read genes
die "Usage: cat gwas.20.BF.sort.ann.txt | ./filterGeneSet.pl genes_file.txt" if $genesFile eq '';

readGmt($genesFile);

# Process STDIN
my $l;
while( $l = <STDIN> ) {
	chomp $l;
	my($bf, $pval, $lltot, $lllr, $llmsa, $id1, $id2, $genes1, $ann1, $genes2, $ann2, $genesShared, $anns1, $anns2 ) = split /\t/, $l;

	# Filter lines
	my($gs1) = hasGenes($genes1);
	my($gs2) = hasGenes($genes2);
	if(($gs1 ne '') && ($gs2 ne '')) { print "$l\t$gs1\t$gs2\n"; }
	
}
