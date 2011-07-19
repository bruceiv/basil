#!/usr/bin/env perl

# copies stdin to stdout, prepending the name given as its first command line 
# argument to a non-empty stream

my $line;
if ($line = <STDIN>) {
	print "@ARGV[0]\n";
	print $line;
	while ($line = <STDIN>) { print $line; }
}
