#!/usr/bin/env perl

# converts a GAP group constructed from its generator list (such as the output 
# of nautyGen.gap) into a list of Basil group generators

my $line = "";

#concatenate input into single line
while (<>) {
	s/^\s*(.+)\s*$/$1/;
	$line .= $_;
}

$line =~ s{^G:=Group\(\[ \(}{};	# remove initial brackets
$line =~ s{\),\s*\(}{\n}g;		# substitute group element delimiters
$line =~ s{\) \]\);;$}{\n};		# remove final brackets

$line =~ s{,}{ }g;				# change element delimeter from comma to space
$line =~ s{\)\(}{ , }g;			# change cycle delimeter from ")(" to " , "

print "symmetry begin\n";
print $line;
print "symmetry end\n";