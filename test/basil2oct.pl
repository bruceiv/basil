#!/usr/bin/env perl

# grabs the matrix and cobasis lines from verbose Basil output, to test with 
# Octave for full rank

while (<>) {
	last if (m/^Matrix/); # skip up to matrix line
}

print "M = [\n";
while (<>) {
	last if (m/^Group/); # end at group line
	
	# strip brackets and print
	s{\[}{}g;
	s{\]}{}g;
	print;
}
print "];\n";

print "C = [\n";
while (<>) {
	if (m/cobases/) {
		# strip surroundings and commas
		s/.*\{//;
		s/\}//;
		s/,//g;
		print;
	}
}
print "];\n";

