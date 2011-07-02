#!/usr/bin/env perl

# grabs the cobasis lines from Symbal output, returning them in sorted order, 
# with the total number of basis orbits printed at the end

my @cobases = ();

while (<>) {
	push(@cobases, $1) if ( m/rec\( cobasis := (\[.*?\]),/ );
	print $_ if ( m/^basis_orbits => \d+,$/ );
}

print join(",\n", sort @cobases), "\n";
