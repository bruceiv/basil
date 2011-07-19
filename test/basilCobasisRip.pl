#!/usr/bin/env perl

# grabs the cobasis lines from Basil output, to compare with Symbal (reformats 
# to GAP output to compare with symbalCobasisRip.pl)

my @cobases = ();
my $norbits = 0;

my $doPrint = 0;
while (<>) {
	if ( $doPrint ) {
		if ( m|\{(.+)\}| ) { push(@cobases, "[".$1."]"); }
		elsif ( m|\}| ) { $doPrint = 0; }
	} elsif ( m/^\tbasis orbits: (\d+)/ ) {
		# start printing at header
		$doPrint = 1;
		$norbits = $1;
	}
}

print "Cobs:=[\n\t";
print join(",\n\t", @cobases);
print "\n];;\n";
print "NOrbits:=${norbits};;\n";
