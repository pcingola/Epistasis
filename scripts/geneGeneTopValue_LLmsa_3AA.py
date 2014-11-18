#!/usr/bin/env python

$lineMax = "";
$max = 0;
for( $count = 0 ; $l = <STDIN> ; ) {
	chomp $l;
	($id1, $id2, $ll, $llnull, $llalt, $seq1, $seq2) = split /\t/, $l;
	
	if( $id1 ne '' ) {
		if( $llalt > $max ) { 
			$max = $llalt;
			$lineMax = $l;
		}

		$count++;
	}
}

print "$count\t$lineMax\n";
