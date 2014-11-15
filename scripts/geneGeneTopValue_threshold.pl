#!/usr/bin/perl

$th = $ARGV[0];
die "Usage: cat gene_gene.txt | ./geneGeneTopValue_threshold.pl LL_null_threshold\n" if $th eq '';

$lineMax = "";
$max = 0;
for( $count = 0 ; $l = <STDIN> ; ) {
	chomp $l;
	($id1, $id2, $ll, $llalt, $llnull, $seq1, $seq2) = split /\t/, $l;
	
	if( $id1 ne '' ) {
		if( $llnull >= $th ) {
			if( $ll > $max ) { 
				$max = $ll;
				$lineMax = $l;
			}

			$countPass++;
		}
		$count++;
	}
}

print "threshold_ll_null: $th\tll_max: $max\tcount: $count\tcount_pass: $countPass\t$lineMax\n";
