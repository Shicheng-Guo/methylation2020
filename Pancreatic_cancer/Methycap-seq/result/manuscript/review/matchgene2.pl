#!/usr/bin/perl
use strict;
use Cwd;
chdir getcwd;

open F,"tumor_related_methylation_abstract.txt";
chomp(my @database=<F>);
close F;


open F1,"suppressionorpromotion.txt";
while(<F1>){
chomp;
print "$_\n" if @database ~~ /\b$_\b/g;
}




