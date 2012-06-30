#!/usr/bin/env perl

# Strips the Gram matrix from a Basil input file

while (<>) {
	last if (m/^gram begin$/); #print lines leading up to gram begin
	print;
}

#skip Gram matrix
while (<>) {
	last if (m/^gram end$/);
}

#print lines after Gram matrix
while (<>) {
	print;
}

