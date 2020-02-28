#!/usr/bin/perl -w

use strict;

my $strand = $ARGV[0];
my $type = $ARGV[1];
my $qual_base = 33;
my %table;

my $list = "/media/2TB_storeA/dbSNP/bspp_snp_positions_12122012";

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

my %chrSizes;
$chrSizes{"chr9"}=141213431;
$chrSizes{"chr5"}=180915260;
$chrSizes{"chr2"}=243199373;
$chrSizes{"chr17"}=81195210;
$chrSizes{"chr14"}=107349540;
$chrSizes{"chr4"}=191154276;
$chrSizes{"chr16"}=90354753;
$chrSizes{"chr21"}=48129895;
$chrSizes{"chrM"}=16571;
$chrSizes{"chr7"}=159138663;
$chrSizes{"chr3"}=198022430;
$chrSizes{"chr18"}=78077248;
$chrSizes{"chr12"}=133851895;
$chrSizes{"chrX"}=155270560;
$chrSizes{"chr13"}=115169878;
$chrSizes{"chr15"}=102531392;
$chrSizes{"chr8"}=146364022;
$chrSizes{"chr22"}=51304566;
$chrSizes{"chr11"}=135006516;
$chrSizes{"chr10"}=135534747;
$chrSizes{"chr20"}=63025520;
$chrSizes{"chr19"}=59128983;
$chrSizes{"chr6"}=171115067;
$chrSizes{"chr1"}=249250621;
$chrSizes{"chrY"}=59373566;

sub main{
	readlist() if($type eq "REF");
	#chr1    11391   t       W       33      33      30      6       ...AAA  ZgccgX
	while(my $line = <STDIN>){
		chomp($line);
		my @f = split /\t/, $line;
		my $c_pos = $chrSizes{$f[0]} - $f[1] + 1;
		next if($f[5] < 20 && $type eq "VAR");
		next if($f[5] > 5 && $type eq "REF");
		next if(!$table{$f[0].":".$c_pos} && $strand eq 'C' && $type eq "REF");
		my %variantStat = maq_pileupFields2variantStat(\@f, 0, $strand, "SNP");
	}
}

sub maq_pileupFields2variantStat(){
        my ($h_fields, $mask_methyl, $str, $context) = @_;
        my @fields = @{$h_fields};
        my $refBase = uc($fields[2]);
        my $readBase = $fields[8];
	my $call = uc($fields[3]);
	$refBase = $rcTable{$refBase} if($str eq 'C');
	#return if($refBase eq 'C' && $str eq 'W' && $call =~ m/[CTY]/);
	#return if($refBase eq 'G' && $str eq 'C' && $call =~ m/[CTY]/);
	#return if($refBase eq 'G' && $str eq 'C' && $call =~ m/[GAR]/);
	#$call = $rcTable{$call} if($str eq 'C');
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
	$fields[1] = $chrSizes{$fields[0]}-$fields[1]+1 if($str eq 'C');

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
sub readlist{
	open(INFILE, "$list") || die("Error opening file\n");
	while(my $line = <INFILE>){
		chomp($line);
		my @f = split /\t/, $line;
		$table{$f[0].":".$f[1]} = 1;
	}
	close(INFILE);
}
main();
