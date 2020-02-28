use strict;
use Cwd;
use List::Util qw( min max );
my $dir=getcwd;
die USAGE() if scalar(@ARGV<1);
my $file=shift @ARGV;
open F,$file || die "cannot file the input $file\n";
my %region;
while(<F>){
next if /Pos|chr/i;
chop;
my ($gene,$pos1,$chr,$pos2,$dis2tss,$type,undef)=split/\t/;
push @{$region{"$gene:$chr"}},$pos2;	
}
foreach my $gene(%region){
my $min = min @($region{$gene});
my $max = max @($region{$gene});
my ($gene,$chr)=split/:/,$gene;
print "chr$chr\t$min\t$max\t$gene\n";
}