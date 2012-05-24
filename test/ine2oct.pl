#!/usr/bin/env perl

# Converts the matrix from a Basil input file into an equivalent Octave file.

while (<>) {
	last if (m/^begin$/); #skip lines leading up to begin
}

#skip dimension line
$line = <>;

print "M = [\n";
while (<>) {
	last if (m/^end$/); #read until end line
	
	print;
}
print "];\n";

