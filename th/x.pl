#!/usr/bin/perl

while( $l = <STDIN> ) {
	
	chomp $l;

	if( $l =~/title=[\"\{](.*)[\"\}]\s*,\s*$/ ) {
		# Parse title
		$title = $1;
		print "TITLE: |$title|\n";
	} else { 
		print "$l\n"; 
	}
}
