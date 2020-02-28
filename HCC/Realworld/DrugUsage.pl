#!/usr/bin/perl
use strict;
chdir "/home/shg047/Dropbox/Realworld";
open F,"DrugUsage.txt" || die "cannot open !";
my %recom;
while(<F>){
chomp;
my ($center,$subject,$visit,$peseq,$drug,$dos)=split/\s+/;
next if $drug=~/!/;
next if $drug eq "";
push @{$recom{"$center-$subject-$visit"}},$drug;
}

foreach my $id(sort keys %recom){
print "$id\t";
foreach my $drug(sort @{$recom{$id}}){
print "$drug;"
}
print "\n";
}
