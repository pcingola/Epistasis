#!/usr/bin/perl

use strict;


#-------------------------------------------------------------------------------
# Simplify annotations
#-------------------------------------------------------------------------------
sub annotate($) {
	my($ann) = @_;

	my($eff, $gene, $impact, $e, $g, @anns, %effs, %genes);
	@anns = split /,/, $ann;
	foreach $ann ( @anns ) {
		# Format example
		# 		NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|cGc/cCc|R433P|714|CAPN1|protein_coding|CODING|ENST00000279247|11|1)
		if( $ann =~ /(.+)\((.*?)\|.*?\|.*?\|.*?\|.*?\|(.*?)\|.*/ ) {
			($eff, $impact, $gene) = ($1, $2, $3);

			if(($impact eq 'MODERATE') || ($impact eq'HIGH')) {
				$effs{$eff} = 1;
				$genes{$gene} = 1;
			}
		}
	}

	my($gstr) = "";
	foreach $g ( sort keys %genes ) { 
		if( $gstr ne '' )	{ $gstr .= ","; }
		$gstr .= $g; 
	}

	my($estr) = "";
	foreach $e ( keys %effs ) { 
		if( $estr ne '' )	{ $estr .= ","; }
		$estr .= $e; 
	}

	return "$gstr\t$estr";
}

#-------------------------------------------------------------------------------
# Share
#-------------------------------------------------------------------------------
sub shared($$) {
	my($ann1, $ann2) = @_;

	my($g, $genes1, $effs1, $genes2, $effs2, %gs, %gs1, %gs2, @gs1, @gs2);
	%gs = ();

	# Parse annotations 1
	($genes1, $effs1) = split /\t/, $ann1;
	@gs1 = split /,/, $genes1; 
	%gs1 = ();
	foreach $g ( @gs1 ) { 
		$gs1{$g} = 1; 
		$gs{$g} = 1; 
	}

	# Parse annotations 2
	($genes2, $effs2) = split /\t/, $ann2;
	@gs2 = split /,/, $genes2; 
	%gs2 = ();
	foreach $g ( @gs2 ) { 
		$gs2{$g} = 2; 
		$gs{$g} = 1; 
	}

	my($gsStr) = '';
	foreach $g ( sort keys %gs ) {
		if( $gs1{$g} ne '' && $gs2{$g} ne '' ) {
			if( $gsStr ne '' )	{ $gsStr .= ','; }
			$gsStr .= $g;
		}
	}

	if( $gsStr eq '' )	{ return "-"; }
	return $gsStr;
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

# Parse input
my($l, $ann1Str, $ann2Str, $share, $i);

my @removeLabels = (3 ,4 ,5 ,6 ,7 ,8 ,9 ,10 ,11 ,12 ,17 ,18);

while( $l = <STDIN> ) {
	chomp $l;
	my @t = split /\t/, $l;

	# Split semicolon and discard label + Remove spaces
	foreach $i (@removeLabels) {
		if( $t[$i] =~/(.+):(.+)/ ) {
			$t[$i] = $2;
			$t[$i] =~ tr/ //ds;
		}
	}

	# Extract genes for MODERATE and HIGH impact variants
	$t[15] = annotate($t[15]);
	$t[16] = annotate($t[16]);

	# Show results
	print join("\t", @t) . "\n";
}
