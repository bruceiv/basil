#!/usr/bin/env perl

# grabs the group generator lines and initial cobasis from Basil output, to 
# analyze symmetries

my @gens = ();
my $initCob = '';

my $readGens = 0;
while (<>) {
	if ( $readGens ) {
		if ( m|(\(.+\))| ) { push(@gens, $1); }
		elsif ( m|\}| ) { $readGens = 0; }
	} elsif ( m|^\tsymmetry generators:| ) {
		# start printing at header
		$readGens = 1;
	} elsif ( m|^\tinitial cobasis: \{(.+)\}| ) {
		$initCob = '['.$1.']';
	}
}

print "InitCob:=${initCob};;\n";
print "Gens:=[\n\t";
print join(",\n\t", @gens);
print "\n];;\n";
