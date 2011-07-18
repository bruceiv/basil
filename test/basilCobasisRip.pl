#!/usr/bin/env perl

# grabs the cobasis lines from Basil output, to compare with Symbal

my $doPrint = 0;
while (<>) {
	# start printing at header
	$doPrint = 1 if m/^\tbasis orbits/;
	print $_ if $doPrint;
	#stop printing at closing bracket
	$doPrint = 0 if $doPrint and ( m/^\t\}/ or m/^\t\{\}/ );
}
