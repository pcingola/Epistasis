#!/usr/bin/perl

#-------------------------------------------------------------------------------
# Extract gene names from annotations
#-------------------------------------------------------------------------------
sub genes($) {
	my($ann) = @_;

	# Split each annotation
	@anns = split /,/, $ann;
	my(%genes);
	foreach $a (@anns) {
		# Extract gene names
		@f = split /\|/, $a;
		$gene = @f[5];
		$genes{$gene} = 1
	}

	# Concatenate all gene names
	my($gs) = "";
	foreach $g ( sort keys %genes ) {
		$gs .= "," if $gs ne '';
		$gs .= $g;
	}

	return $gs;
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

while( $l = <STDIN> ) {
	chomp $l;

	@t = split /\t/, $l;
	$anni = @t[6];
	$annj = @t[7];
	$genei = genes($anni);
	$genej = genes($annj);
	print "$genei\t$genej\t$l\n";
	
}
