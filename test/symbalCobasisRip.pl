#!/usr/bin/env perl

# grabs the cobasis lines from Symbal output, returning them in sorted order, 
# with the total number of basis orbits printed at the end. This is also 
# usable GAP code, for ease of automated comparison

my @cobases = ();
my $linebuf = "";
my $line = "";
my $norbits = 0;

while (<>) {
	# get all cobases
	# push(@cobases, "{".$1."}") if ( m/rec\( cobasis := \[(.*?)\],/ );
	if ( m|\s+cobasis\s+| ) {
		$linebuf = $_;
		while ( ! ( $linebuf =~ m|\s+cobasis := \[(.*?)\],| ) ) {
			#add another line to the input
			$line = <>; $line =~ s|^\s*||;
			$linebuf =~ s|\s*$||;
			$linebuf = $linebuf.$line;
		}
		$linebuf =~ m|\s+cobasis := \[(.*?)\],|;
		push(@cobases, "[".$1."]");
	}
	# print basis orbits header
# 	print "\tbasis orbits: ", $1, "\n" if ( m/^basis_orbits => (\d+),$/ );
	$norbits = $1 if ( m/^basis_orbits => (\d+),$/ );
}

# sort cobases (uses a lexicographical sort, but doesn't have to account for 
# different lengths, or identical elements)
@cobases = sort {
	my $cv = 0;
	# reset match position for variables
	pos($a) = -1; pos($b) = -1;
	while ( $cv == 0 ) {
		#match the next number
		$a =~ m/\d+/g; my $av = $&;
		$b =~ m/\d+/g; my $bv = $&;
		#and compare numerically
		$cv = $av <=> $bv;
	}
	$cv;
} @cobases;

print "Cobs:=[\n\t";
print join(",\n\t", @cobases);
print "\n];;\n";
print "NOrbits:=${norbits};;\n";
