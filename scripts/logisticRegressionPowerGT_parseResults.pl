#!/usr/bin/perl

$pvalTh = 5E-14;

print "perc\tn\taf1\taf2\taf12\tbeta3\tcount\ttot\n";
$tot = 0;
while( $l = <STDIN> ) {
	chomp $l;

	if( $l =~ /^#/ ) {
		if( $l =~ /^#---/ ) {

			# Show previous stats
			if( $tot > 0 ) {
				$perc = 100 * $count / $tot;
				print "$perc\t$n\t$af1\t$af2\t$af12\t$beta3\t$count\t$tot\n";
			}
			$count = $tot = 0;
		}
	} else {
		# Parse line
		$l =~ tr/:/\t/;
		@t = split /\t/, $l;

		if( $#t > 12 ) {
			$n = $t[2];
			$af1 = $t[4];
			$af2 = $t[6];
			$af12 = $t[8];
			$beta3 = $t[10];
			$pval = $t[12];
        
			if( $pval < $pvalTh ) { $count++; }
			$tot++;
		}
	}
}

# Show last  stats
if( $tot > 0 ) {
	$perc = 100 * $count / $tot;
	print "$perc\t$n\t$af1\t$af2\t$af12\t$beta3\t$count\t$tot\n";
}
