#!/usr/bin/perl

#---
# Initialize AA map
#---

$n = 0;
$numByAa{'A'} = $n; $aaByNum[$n] = 'A'; $n++;
$numByAa{'C'} = $n; $aaByNum[$n] = 'C'; $n++;
$numByAa{'D'} = $n; $aaByNum[$n] = 'D'; $n++;
$numByAa{'E'} = $n; $aaByNum[$n] = 'E'; $n++;
$numByAa{'F'} = $n; $aaByNum[$n] = 'F'; $n++;
$numByAa{'G'} = $n; $aaByNum[$n] = 'G'; $n++;
$numByAa{'H'} = $n; $aaByNum[$n] = 'H'; $n++;
$numByAa{'I'} = $n; $aaByNum[$n] = 'I'; $n++;
$numByAa{'K'} = $n; $aaByNum[$n] = 'K'; $n++;
$numByAa{'L'} = $n; $aaByNum[$n] = 'L'; $n++;
$numByAa{'M'} = $n; $aaByNum[$n] = 'M'; $n++;
$numByAa{'N'} = $n; $aaByNum[$n] = 'N'; $n++;
$numByAa{'P'} = $n; $aaByNum[$n] = 'P'; $n++;
$numByAa{'Q'} = $n; $aaByNum[$n] = 'Q'; $n++;
$numByAa{'R'} = $n; $aaByNum[$n] = 'R'; $n++;
$numByAa{'S'} = $n; $aaByNum[$n] = 'S'; $n++;
$numByAa{'T'} = $n; $aaByNum[$n] = 'T'; $n++;
$numByAa{'V'} = $n; $aaByNum[$n] = 'V'; $n++;
$numByAa{'W'} = $n; $aaByNum[$n] = 'W'; $n++;
$numByAa{'Y'} = $n; $aaByNum[$n] = 'Y'; $n++;

#---
# Load matrix
#---
while( $l = <STDIN> ) { 
	chomp $l;
	($aa1, $aa2, $val) = split /\t/, $l;
	$d[ $numByAa{$aa1} ][ $numByAa{$aa2} ] = $val;
}

#---
# Show as matrix
#---

# Title
for( $i=0 ; $i < $n ; $i++ ) {
	print "$aaByNum[$i]";
	print "\t" if( $i < ($n-1) );
}
print "\n";

# Entries
for( $i=0 ; $i < $n ; $i++ ) {
	print "$aaByNum[$i]";
	for( $j=0 ; $j < $n ; $j++ ) {
		if( $d[$i][$j] )	{ print "\t$d[$i][$j]"; }
		else				{ print "\t0"; }
	}
	print "\n";
}
