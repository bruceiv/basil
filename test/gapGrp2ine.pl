#!/usr/bin/env perl

my $line = "";

#concatenate input into single line
while (<>) {
	s/^\s*(.+)\s*$/$1/;
	$line .= $_;
}

$line =~ s{^\[ \(}{};		# remove initial list bracket
$line =~ s{\), \(}{\n}g;	# substitute group element delimiters for newlines
$line =~ s{\) \]$}{\n};		# remove final list bracket

$line =~ s{,}{ }g;			# change element delimeter from comma to space
$line =~ s{\)\(}{ , }g;		# change cycle delimeter from ")(" to " , "

print "symmetry begin\n";
print $line;
print "symmetry end\n";