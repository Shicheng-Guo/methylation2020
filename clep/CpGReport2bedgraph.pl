#!/usr/bin/perl
use strict;
use Cwd;
my $file=shift @ARGV;
open F,"gunzip -c $file|" || die "cannot open the $file\n";
open OUT,">$file.bdg";
while(<F>){
my $c=$_;
my $g=<F>;
# print "$c";
# print "$g";
my @c=split/\s+/,$c;
my @g=split/\s+/,$g;

if($c[2]=="+" and $g[2]=="-" and $g[1]=$c[1]+1){
my $mf=sprintf "%.2f",($c[3]+$g[3])/($c[3]+$g[3]+$c[4]+$g[4]+0.1);
print OUT "$c[0]\t$c[1]\t$g[1]\t$mf\n";
}
}
close OUT;

