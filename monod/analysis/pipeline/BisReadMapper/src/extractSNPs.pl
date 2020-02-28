#!/usr/bin/perl -w

use strict;

my $strand = $ARGV[0];
my $type = $ARGV[1];
my $qual_base = $ARGV[2];
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
	
	#chr1    11391   t       W       33      33      30      6       ...AAA  ZgccgX
	while(my $line = <STDIN>){
		chomp($line);
		my @f = split /\t/, $line;
		next if($f[5] < 20 && $type eq "VAR");
		next if($f[5] > 5 && $type eq "REF");
		my %variantStat = maq_pileupFields2variantStat(\@f, 0, $strand, "SNP");
	}
}

sub maq_pileupFields2variantStat(){
        my ($h_fields, $mask_methyl, $str, $context) = @_;
        my @fields = @{$h_fields};
        my $refBase = uc($fields[2]);
        my $readBase = $fields[8];
	my $call = uc($fields[3]);
	$call = $rcTable{$call} if($str eq 'C');
	return if(!$call);
        my $readQual = $fields[9];
        $readBase =~ s/\$//g;
        $readBase =~ s/\^.//g;
        $readBase =~ s/F//g;
        $readBase =~ s/-[0-9]+[ACGTNacgtn]+/D/g;
        $readBase =~ s/\+[0-9]+[ACGTNacgtn]+/I/g;
        #$readBase =~ s/[\.\,]/$refBase/g;
        my $totalCounts=0;
        my %variantStat;
        while(my $base = chop($readBase)){
                my $baseQual = chop($readQual);
		#print $base, ":", ord($baseQual)-$qual_base, " ";
                next if(ord($baseQual)-$qual_base < 5); #ignore low-quality bases (phred score < 5)
		next if($base !~ /[ATGC\.]/ && $str eq 'W'); #ignore not uppercase letters and dots if counting Watson
		next if($base !~ /[atgc\,]/ && $str eq 'C'); #ignore not lowercase letters and commas if counting Crick
		$base =~ s/[\.\,]/$refBase/;
		$base = uc($base); # make all uppercase
		$base = $rcTable{$base} if($str eq 'C');
		#print $base, " ";
                $base =~ s/C/T/g if($mask_methyl); # mask methylation
                $variantStat{'counts'}->{$base}++;
                $totalCounts++;
        }
	#print "\n";
        $variantStat{'refBase'}=$refBase;
        $variantStat{'depth'}= $totalCounts;
        $variantStat{'snpQual'}= $fields[5];
        $variantStat{'call'}= $call;

	return if($variantStat{'depth'} == 0);
	print $fields[0], "\t", $fields[1], "\t", 
		$variantStat{'refBase'}, "\t", $str, "\t", 
		$variantStat{'call'}, "\t", $variantStat{'snpQual'}, "\t",  $variantStat{'depth'};
	foreach my $base (keys(%{$variantStat{'counts'}})){
		my $cnt = $variantStat{'counts'}->{$base};
		$base = $rcTable{$base} if($str eq 'C');
		print "\t$base\t$cnt";
        }
        print "\n";
        return %variantStat;
}

main();
