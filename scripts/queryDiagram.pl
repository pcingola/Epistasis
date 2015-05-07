#!/usr/bin/perl

#-------------------------------------------------------------------------------
# Find closes marker and show data
#-------------------------------------------------------------------------------
sub find($$) {
	($chr, $pos) = @_;
	for( $i=0 ; $i < 1000000 ; $i++ ) {
		$p = $pos + $i;
		$key = "$chr:$p";

		if( $pvals{$key} ne '' ) {
			print "$chr:$pos\t$i\t$key\t'$pvals{$key}'\t'$ors{$key}'\n";
			return;
		}

		$p = $pos - $i;
		$key = "$chr:$p";
		if( $pvals{$key} ne '' ) {
			print "$chr:$pos\t$i\t$key\t'$pvals{$key}'\t'$ors{$key}'\n";
			return;
		}
	}

	print "$chr:$pos\tNOTHING FOUND!\n";
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

# Debug
$debug = 0;

# Read file
print STDERR "Reading:  DIAGRAMv3.2012DEC17.txt\n";
open D,"DIAGRAMv3.2012DEC17.txt" || die "Cannot open file\n";
while( $l = <D> ) {
	chomp $l;
	($snpId,$chr,$pos,$alt,$ref,$pval,$or,$or_95l,$or_95u) = split /\t/, $l;
	#$key = "$chr:$pos\_$ref/$alt";
	$key = "$chr:$pos";

	$pvals{$key} .= "$pval\t";
	$ors{$key} .= "$or [$or_95l , $or_95u]\t";
}
close D;

# Read file
print STDERR "Reading: DIAGRAMv3.2013MAY07.metabochip.txt\n";
open D,"DIAGRAMv3.2013MAY07.metabochip.txt" || die "Cannot open file\n";
while( $l = <D> ) {
	chomp $l;
	($chr,$pos,$snpId,$alt,$ref,$or,$or_95l,$or_95u,$pval) = split /\t/, $l;
	#$key = "$chr:$pos\_$ref/$alt";
	$key = "$chr:$pos";

	$pvals{$key} .= "$pval\t";
	$ors{$key} .= "$or [$or_95l , $or_95u]\t";
}
close D;

# Read file
print STDERR "Reading: DIAGRAMv3.2014OCT16.TransEthnic_T2D_GWAS.MegaMeta.txt\n";
open D,"DIAGRAMv3.2014OCT16.TransEthnic_T2D_GWAS.MegaMeta.txt" || die "Cannot open file\n";
while( $l = <D> ) {
	chomp $l;
	($snpId,$chr,$pos,$alt,$ref,$or,$or_95l,$or_95u,$pval) = split /\t/, $l;
	#$key = "$chr:$pos\_$ref/$alt";
	$key = "$chr:$pos";

	$pvals{$key} .= "$pval\t";
	$ors{$key} .= "$or [$or_95l , $or_95u]\t";
}
close D;

# Queries
print STDERR "Readinf queries form STDIN\n";
while( $l = <STDIN> ) {
	chomp $l;
	($chr,$pos,$ref,$alt) = split /\t/, $l;
	#$key = "$chr:$pos\_$ref/$alt";
	$key = "$chr:$pos";
	find($chr,$pos)
}

