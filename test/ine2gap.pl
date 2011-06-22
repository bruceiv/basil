#!/usr/bin/env perl

my $vrep = "false";		#boolean V-representation flag
my @mrows = ();			#matrix rows (GAP format)
my @grows = ();			#group rows (GAP format)
my $line;				#temp variable

while (<>) {
	last if (m/^begin$/); #skip lines leading up to begin
	
	#handle V-representation option
	$vrep = "true" if (m/^V-representation$/);
}

#skip dimension line
$line = <>;

while (<>) {
	last if (m/^end$/); #read until end line
	
	# get whitespace separated coordinates
	# NOTE that this assumes one vector per line, tighter than assumed by LRS
 	s/^\s+//; #trim leading whitespace from line
 	my @coords = split(/\s+/);
	
	#convert to GAP list representation
	$line = join(',', @coords);
	$line = '['.$line.']';
	
	push(@mrows, $line);
}

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

print '
RequirePackage("require_file");
AddToLoadPath("~/research/symbal/src");

ReadFile("dfs.gap");

';

print "M:=[";
print join(",\n", @mrows);
print "];;\n\n";

print "G_gen:=[";
print join(",\n", @grows);
print "];;\n";

print"G:=Group(G_gen);;\n";

print '
results:=dfs(M,G,[rec(',"VRepresentation:=$vrep",')]);
FormatResults(results,"dimension");

';