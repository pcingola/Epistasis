#!/usr/bin/perl

$num = 100;
$height = 2543;
$alpha = $height / $num;
$beta = 15;

print STDERR "alpha = $alpha\n";
print "<img src=tree.png>\n";
print "<svg height=2543 width=200>\n";

print "\t<text x=0 y=1>---------</text>\n";
for( $i=0 ; $l = <STDIN>; $i++ ) {
	chomp $l;

	if( $l !~ /-/ ) {
		$y = $alpha * $i + $beta;
		print "\t<text x=0 y=$y>$l</text>\n";
	}
}

print "</svg>\n";
