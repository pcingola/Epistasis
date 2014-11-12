#!/usr/bin/perl

$debug = 0;

#---
# Read clinvar
#---
open CLINVAR, "../clinvar-latest.eff.vcf" or die "Cannot open file\n";
while( $l = <CLINVAR> ) {
	chomp $l;
	($chr, $pos, $id, $ref, $alt, $q, $filter, $info) = split /\t/, $l;

	if( $info =~/CLNSIG=(.*?);/ ) {
		$clnsig = $1;
		if( $clnsig =~ /(\d+)[,|].*/ )	{ 
			print "$clnsig\t=>\t$1\n" if $debug;
			$clnsig = $1; 
		}

		$clnsig{"$chr:$pos"} = $clnsig;
		print "$chr\t$pos\t$clnsig\n" if $debug;
	}
}
close CLINVAR;

# Read STDIN
print "chr\tpos\tref\talt\tll\tlen\tclnsig\n";
while( $l = <STDIN> ) {
	chomp $l;
	($chr, $pos, $ref, $alt, $ll, $len) = split /\t/, $l;

	$pos += 1;
	print "$chr\t$pos\t$ref\t$alt\t$ll\t$len\t" . $clnsig{"$chr:$pos"} . "\n";
}
