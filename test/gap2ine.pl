#!/usr/bin/env perl

use Text::Balanced qw(extract_bracketed);

undef $/;

$in=<>;

$in=~s|\s+||g;

$in=~ s|^\s*(\w+):=||;
$in=~s|(-?[0-9]+\s*/\s*[0-9]*)|'\1'|g;

$aref=eval($in);

my $rows=scalar(@{$aref});
my $cols=scalar(@{$aref->[1]});

print "begin\n$rows $cols rational\n";
for (my $i=0; $i<$rows; $i++){
  for (my $j=0; $j<$cols; $j++){
    print $aref->[$i]->[$j], " ";
  }
  print "\n";
}
print "end\n";
print "geometric\n";
