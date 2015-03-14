#!/usr/bin/perl

#-------------------------------------------------------------------------------
# Find interactions appearing in biogrid
#
#															Pablo Cingolani 2015
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

# Parse command line parameters
$bioGridFile = $ARGV[0];

#---
# Read entries form biogrid
#---

print STDERR "Reading $bioGridFile\n";
open BIOGRID, $bioGridFile or die "Cannot open biogrid file '$bioGridFile'\n";
while( $l = <BIOGRID> ) {
	chomp $l;
	($g1, $g2) = split /\t/, $l;

	# Add interaction (both ways)
	$biogrid{"$g1\t$g2"} = 1;
	$biogrid{"$g2\t$g1"} = 1;
}
print STDERR "done.\n";

#---
# Read STDIN (results form GWAS analysis)
#---

while( $l = <STDIN> ) {
	chomp $l;
	@t = split /\t/, $l;
	
	($g1, $g2) = ($t[15], $t[17]);
	$key = "$g1\t$g2";

	
	if( $biogrid{$key} )	{ print "$l\tBIOGRID\n"; }
	else					{ print "$l\t\n"; }
}

