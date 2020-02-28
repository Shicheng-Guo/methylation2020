#!/usr/bin/perl -w
#Usage: ./extractMethyl.pl [position list] [quality PHRED base] < [vcf file]

use strict;

my $minPHRED = 5;
my $cpg_list = $ARGV[0];
my $qual_base = $ARGV[1];
my %cpgTable;

my %rcTable;
$rcTable{'A'}='T';
$rcTable{'T'}='A';
$rcTable{'G'}='C';
$rcTable{'C'}='G';
$rcTable{'N'}='N';
$rcTable{'R'}='Y';
$rcTable{'Y'}='R';
$rcTable{'M'}='K';
$rcTable{'K'}='M';
$rcTable{'S'}='S';
$rcTable{'W'}='W';

sub main{
	#chr10:W 60025   CG
	open(CPG_LIST, "$cpg_list") || die("Error opening $cpg_list\n");
	while(my $line = <CPG_LIST>){
		chomp($line);
		my @f = split /\t/, $line;
		my ($chr, $str) = split ":", $f[0];
		my $index = $chr.":".$f[1];
		$cpgTable{$index}=$f[2].":".$str;
	}
	close(CPG_LIST);

	#chr1    11312   g       3       aaa     dO^
	while(my $line = <STDIN>){
		chomp($line);
		my @f = split /\t/, $line;
		my $index = $f[0]. ":". $f[1];
		if($cpgTable{$index}){
			my ($context, $str) = split ":", $cpgTable{$index};
			print "Unexpected base at $index\n" if($f[2] !~ m/c/i && $str eq "W");
			print "Unexpected base at $index\n" if($f[2] !~ m/g/i && $str eq "C");
			my %variantStat = pileupFields2variantStat(\@f, 0, $str, $context);
		}else{
			$f[1]--;
			my $indexC = $f[0]. ":" . $f[1];
			if($cpgTable{$indexC}){
				#this is on CG with C on the Crick strand, only applies to CG positions
				my ($context, $str) = split ":", $cpgTable{$indexC};
				next if($context ne "CG");
				print "Unexpected base at $index\n" if($f[2] !~ m/g/i);
				my %variantStat = pileupFields2variantStat(\@f, 0, "C", $context);
			}
		}	
	}
}

sub pileupFields2variantStat(){
        my ($h_fields, $mask_methyl, $str, $context) = @_;
        my @fields = @{$h_fields};
        my $refBase = uc($fields[2]);
        my $readBase = $fields[4];
        my $readQual = $fields[5];
        my %variantStat;
	if(!$readBase || !$readQual){ return %variantStat;}
        $readBase =~ s/\$//g;
        $readBase =~ s/\^.//g;
        $readBase =~ s/F//g;
	$readBase =~ s/S//g;
        my $totalCounts=0;
	#print $readBase, "\n";
	#print $readQual, "\n"; 
	my @quals = split "", $readQual;
	my $j = 0;
	my $i = 0;
	while($i < length($readBase)){
		if(substr($readBase, $i, 1) =~ /^[+-]/){
			#skip indels
			substr($readBase,$i+1) =~ /^([0-9]+)/;
			defined($1) || die("Could not understand indel at pos $i: $readBase\n");
			my $run_len = length($1)+$1+1;
			$i+=$run_len;
		}else{
			my $base = substr($readBase,$i,1);
			my $baseQual = $quals[$j];
			defined($baseQual) || die("Could not get quality for $j at $i: \n $readBase\n $readQual\n");
			$j++;
			$i++;
			next if(ord($baseQual)-$qual_base < $minPHRED); #ignore low-quality bases (phred score < 5)
			next if($base !~ /[ATGC\.]/ && $str eq 'W'); #ignore not uppercase letters and dots if counting Watson
			next if($base !~ /[atgc\,]/ && $str eq 'C'); #ignore not lowercase letters and commas if counting Crick
			$base =~ s/[\.\,]/$refBase/;
			$base = uc($base); # make all uppercase
			$base =~ s/C/T/g if($mask_methyl); # mask methylation
			$variantStat{'counts'}->{$base}++;
			$totalCounts++;
		}
	}

	$variantStat{'refBase'}=$refBase;
        $variantStat{'depth'}= $totalCounts;
	return if($variantStat{'depth'} == 0);
	print $fields[0], "\t", $fields[1], "\t", $str, "\t", $variantStat{'depth'}, "\t", $context;
	foreach my $base (keys(%{$variantStat{'counts'}})){
		my $cnt = $variantStat{'counts'}->{$base};
		$base = $rcTable{$base} if($str eq 'C');
		print "\t$base\t$cnt";
        }
        print "\n";
        return %variantStat;
}

main();
