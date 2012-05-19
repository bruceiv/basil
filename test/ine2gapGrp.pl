#!/usr/bin/env perl

# Gets the GAP group from a Basil input file

my @grows = ();			#group rows (GAP format)
my $line;				#temp variable

#skip dimension line
$line = <>;

while (<>) {
	#skip option lines until beginning of group specification
	last if (m/^symmetry begin$/);
}

while (<>) {
	last if (m/^symmetry end$/); #read until end of group specification
	
	#split cycles on ','
	my @cycles = split(',');
	$line = '';
	
	for my $cycle (@cycles) {
		# get whitespace separated elements
		$cycle =~ s/^\s+//; #trim leading whitespace
		my @elems = split(/\s+/, $cycle);
		
		my $perm = join(',', @elems);
		$perm = '('.$perm.')';
		
		$line .= $perm;
	}
	
	push(@grows, $line);
}

print "G_gen:=[";
print join(",\n", @grows);
print "];;\n";

print"G:=Group(G_gen);;\n";
