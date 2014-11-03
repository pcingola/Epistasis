#!/usr/bin/perl

use strict;

#-------------------------------------------------------------------------------
# Simplify annotations
#-------------------------------------------------------------------------------
sub annotate($) {
	my($ann) = @_;

	my($eff, $gene, $e, $g, @anns, %effs, %genes);
	@anns = split /,/, $ann;
	foreach $ann ( @anns ) {
		# Format example
		# 		NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|cGc/cCc|R433P|714|CAPN1|protein_coding|CODING|ENST00000279247|11|1)
		if( $ann =~ /(.+)\(.*\|.*\|.*\|.*\|.*\|(.*)\|.*\|.*\|.*\|.*\|.*\)/ ) {
			($eff, $gene) = ($1, $2);
			$effs{$eff} = 1;
			$genes{$gene} = 1;
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

my($l, $ann1Str, $ann2Str, $share, $bf, $pva, $ll, $lllr, $llmsa, $id1, $id2, $ann1, $ann2);
while( $l = <STDIN> ) {
	chomp $l;
	#($bf, $pva, $ll, $lllr, $llmsa, $id1, $id2, $ann1, $ann2) = split /\t/, $l;
	($pva, $ll, $lllr, $llmsa, $id1, $id2, $ann1, $ann2) = split /\t/, $l;

	$ann1Str = annotate($ann1);
	$ann2Str = annotate($ann2);

	$share = shared($ann1Str, $ann2Str);

	#print "$bf\t$pva\t$ll\t$lllr\t$llmsa\t$id1\t$id2\t$ann1Str\t$ann2Str\t$share\t$ann1\t$ann2\n";
	print "$pva\t$ll\t$lllr\t$llmsa\t$id1\t$id2\t$ann1Str\t$ann2Str\t$share\t$ann1\t$ann2\n";

}
