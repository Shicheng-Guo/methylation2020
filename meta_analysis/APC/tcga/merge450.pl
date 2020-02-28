#!/usr/bin/perl
use strict;
use Cwd;
chdir "/home/gsc/Dropbox/Project/APCmeta/tcga";
open OUT1, ">fileList.txt";
open OUT2, ">geneList.txt";
#my @cancer=qw(BRCA CESC COAD DLBC ESCA GBM HNSC KICH KIRC KIRP LAML LGG LIHC LUAD LUSC PAAD PRAD READ SARC SKCM STAD THCA UCEC);
my @cancer=qw(LUAD LUSC);
foreach my $cancer(@cancer){
open OUT3, ">$cancer.450+27k.Integrate.txt";
my @file=glob("*edu_$cancer*TCGA*.txt");
my %data;
my $i;
my %version;


foreach my $file(@file){
my @tmp=split/\.|-/,$file;
my $sample=join("-",@tmp[7..10]);
my $version=@tmp[4];
push (@{$version{$sample}},$version);
print "$sample\n";
$i++;
print "$i\n";
open F,$file;
while(<F>){
chomp;
next if /Composite/;
next if /Hybridization/;
my ($cg,$value,$geneid,$chr,$pos)=split/[\||\t]/;
#print "$cg\t$value\n";
$value = sprintf("%.4f",$value) if ($value ne "NA");
$data{$sample}{$cg}{$version}=$value;
#print "$value\n";
}
}

open F2,"share6cpg.txt";
chomp(my @sharecpg=<F2>);
close F2,

my @file= keys %data;
print "$file[0]\n";
my @gene= keys %{$data{$file[0]}};
# print header(file name)
print OUT3 "\tsample\tversion\t";
foreach my $gene(@sharecpg){
print OUT3 "$gene\t";
}
print OUT3 "\n";
foreach my $file(@file){
print "@{$version{$file}}";
my @sortversion=sort {$a<=>$b} @{$version{$file}};
my $maxversion=$sortversion[-1];
print "$maxversion\n";
print OUT3 "$file\t";
print OUT3 "$file\t$maxversion\t";
foreach my $gene(@sharecpg){
print OUT3 "$data{$file}{$gene}{$maxversion}\t";
}
print OUT3 "\n";
}
}
