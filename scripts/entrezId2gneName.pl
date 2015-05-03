#!/usr/bin/perl

#-------------------------------------------------------------------------------
#
# Convert ENTREZ IDs to Gene names
#
#																Pablo Cingolani
#-------------------------------------------------------------------------------

# Parse command line options
$g2e = $ARGV[0];
die "Usage: cat file.txt | entrezId2gneName.pl entrezId2gneName.pl\n" if( $g2e eq '' );

#---
# Read mapping file
#---
print STDERR "Reading IDs from file '$g2e'\n";
open G2E, $g2e;
for( $count=0 ; $l = <G2E> ; $count++ ) {
	chomp $l;
	($gene, $id) = split /\t/, $l;
	$e2g{$id} = $gene;
}
close G2E;
print STDERR "Done: $count lines added\n";

#---
# Read STDIN
#---
while( $l = <STDIN> ) {
	chomp $l;
	(@ids) = split /,/, $l;

	# Convert to name and add to list
	# Note: Each line can have more than two entries
	@genes = ();
	for( $i=0 ; $i <= $#ids ; $i++ ) {
		$gene = $e2g{ $ids[$i] };
		if( $gene ne '' ) {
			push @genes, $gene;
		}
	}

	# Show all possible combinations
	if( $#genes >= 1 ) {
		for( $i=0 ; $i <= $#genes ; $i++ ) {
			for( $j=$i+1 ; $j <= $#genes ; $j++ ) {
				print "$genes[$i]\t$genes[$j]\n";
			}
		}
	}
}
