
#!/usr/bin/perl
use strict;
use Cwd;

my $input= shift @ARGV;
my $referGeneFile="/gpfs/home/guosa/hpc/db/hg19/GeneSymbol2AllAliasName.txt";
my $dir=getcwd;

open F,$referGeneFile;
my %dbgene;
while(<F>){
        chomp;
        my @tmp=split /\s+/;
        my $symbol=$tmp[0];
        $dbgene{$symbol}=$symbol;
}
close F;

open F,$input;
while(<F>){
        chomp;
        my $line=$_;
        foreach my $symbol(keys %dbgene){
          if($line =~ m/\b$symbol\b/i){
          print "$symbol\t$_\n";
          }
        }
}
