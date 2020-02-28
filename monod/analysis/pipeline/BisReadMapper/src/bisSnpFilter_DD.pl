#!/usr/bin/perl -w
# bisSnpFilter.pl: a perl script to combine raw SNP calls from both strands of DNA 
#                   and report confident SNP calls from bisulfite sequencing data. 
# USAGE: bisSnpFilter_DD.pl raw_snp_file dbSnp_file > filtered_snp_file
#

use strict;

my $script_dir = $ARGV[0];
#rules table if reference SNP available
my $RsRules = $script_dir . "/rs_rules";
#rules table if novel SNP
my $DsRules = $script_dir . "/ds_rules";

my %rcTable;
$rcTable{"A"}="T";
$rcTable{"T"}="A";
$rcTable{"G"}="C";
$rcTable{"C"}="G";
$rcTable{"N"}="N";
$rcTable{"R"}="Y";
$rcTable{"Y"}="R";
$rcTable{"M"}="K";
$rcTable{"K"}="M";
$rcTable{"S"}="S";
$rcTable{"W"}="W";
my %three2one;
$three2one{"A/G"}="R";
$three2one{"G/A"}="R";
$three2one{"C/T"}="Y";
$three2one{"T/C"}="Y";
$three2one{"A/C"}="M";
$three2one{"C/A"}="M";
$three2one{"G/T"}="K";
$three2one{"T/G"}="K";
$three2one{"C/G"}="S";
$three2one{"G/C"}="S";
$three2one{"A/T"}="W";
$three2one{"T/A"}="W";
my %hWRules;
$hWRules{"R"}->{"A"} = "R";
$hWRules{"R"}->{"G"} = "R";
$hWRules{"M"}->{"A"} = "M";
$hWRules{"M"}->{"C"} = "M";
$hWRules{"K"}->{"G"} = "K";
$hWRules{"K"}->{"T"} = "K";
$hWRules{"S"}->{"G"} = "S";
$hWRules{"S"}->{"C"} = "S";
$hWRules{"W"}->{"A"} = "W";
$hWRules{"W"}->{"T"} = "W";
my %hCRules;
$hCRules{"Y"}->{"C"} = "Y";
$hCRules{"Y"}->{"T"} = "Y";
$hCRules{"M"}->{"C"} = "M";
$hCRules{"M"}->{"A"} = "M";
$hCRules{"K"}->{"T"} = "K";
$hCRules{"K"}->{"G"} = "K";
$hCRules{"S"}->{"C"} = "S";
$hCRules{"S"}->{"G"} = "S";
$hCRules{"W"}->{"T"} = "W";
$hCRules{"W"}->{"A"} = "W";
my %callRuleW;
my %callRuleC;
my %callRuleDs;
my %candidate_snp_info;

sub main(){
	if(!$ARGV[0]){
		print "Usage: \n";
		exit 0;
	}
	read_snp_list($ARGV[1]);
	load_dbsnp($ARGV[2]);
	call_genotype();
	report_SNPs();
}
sub report_SNPs(){
	#print "SNP position\tSNP call\tSNP qual\tdbSNP\tRefAlleles\tSNP call(fwd)\tAllele count(fwd)\tSNP call(rev)\tAllele count(rev)\n";
	foreach my $index(keys(%candidate_snp_info)){
		next if(!$candidate_snp_info{$index}->{"call"});
		print "$index\t", $candidate_snp_info{$index}->{"call"}, "\t", $candidate_snp_info{$index}->{"qual"}, "\t";
		if($candidate_snp_info{$index}->{"rs"}){
			print $candidate_snp_info{$index}->{"rs"}, "\t", $candidate_snp_info{$index}->{"allele"}, "\t";
		}else{
			print "-\t-\t";
		}
		if($candidate_snp_info{$index}->{"W"}->{"call"}){
			print $candidate_snp_info{$index}->{"W"}->{"call"}, "\t", $candidate_snp_info{$index}->{"W"}->{"allele_count"}, "\t";
		}else{
			print "-\t-\t";
		}
		if($candidate_snp_info{$index}->{"C"}->{"call"}){
			print $candidate_snp_info{$index}->{"C"}->{"call"}, "\t", $candidate_snp_info{$index}->{"C"}->{"allele_count"}, "\t";
		}else{
			print "-\t-\t";
		}
		print "\n";
	}
}

sub read_snp_list(){
	my $filename = shift;
	open(INFILE, "$filename") || die("Error in opening file $filename\n");
	while(my $line = <INFILE>){
		chop($line);
		my @fields = split(/\t/, $line);
		my $key = $fields[0] . ":" . $fields[1];
		my $last_pos = scalar(@fields)-1;
		$candidate_snp_info{$key}->{$fields[3]}->{"refBase"} = $fields[2];
		$candidate_snp_info{$key}->{$fields[3]}->{"call"} = $fields[4];
		$candidate_snp_info{$key}->{$fields[3]}->{"call"} = $rcTable{$fields[4]} if($fields[3] eq "C");
		$candidate_snp_info{$key}->{$fields[3]}->{"qual"} = $fields[5];
		$candidate_snp_info{$key}->{$fields[3]}->{"depth"} = $fields[6];
		$candidate_snp_info{$key}->{"totaldepth"} += $fields[6];
		$candidate_snp_info{$key}->{$fields[3]}->{"allele_count"} = join(",", @fields[7..$last_pos]) if(scalar(@fields)>6);
	}
	close(INFILE);
}

sub call_genotype(){
	open(INFILE, "$RsRules") || die("Error opening file\n");
	while(my $line = <INFILE>){
		chomp($line);
		next if($line =~ m/Call/);
		my @g = split /\t/, $line;
		$callRuleW{$g[0]}->{$g[2]} = $g[4];
		$callRuleC{$g[0]}->{$g[2]} = $g[6];
	}
	close(INFILE);
	open(INFILE, "$DsRules") || die("Error opening file\n");
	while(my $line = <INFILE>){
		chomp($line);
		next if($line =~ m/Call/);
		my @g = split /\t/, $line;
		$callRuleDs{$g[0]}->{$g[2]} = $g[4];
	}
	close(INFILE);
	foreach my $index(keys(%candidate_snp_info)){
		my ($genotype1, $genotype2, $call, $qual);
		my $totaldepth = $candidate_snp_info{$index}->{"totaldepth"};
		next if($totaldepth<8);	
		my $allele = $candidate_snp_info{$index}->{"allele"};
		my $call_w = $candidate_snp_info{$index}->{"W"}->{"call"};
		my $call_c = $candidate_snp_info{$index}->{"C"}->{"call"};
		my $count_w = $candidate_snp_info{$index}->{"W"}->{"allele_count"};
		my $count_c = $candidate_snp_info{$index}->{"C"}->{"allele_count"};
		if($call_w && $call_w =~ m/[AGCT]/){
			$call_w = recalibrate('W', $call_w, $allele, $count_w);
			$candidate_snp_info{$index}->{"W"}->{"call"} = $call_w;
		}
		if($call_c && $call_c =~ m/[AGCT]/){
			$call_c = recalibrate('C', $call_c, $allele, $count_c);
			$candidate_snp_info{$index}->{"C"}->{"call"} = $call_c;
		}
		if($call_w && $call_c){
			$qual = $candidate_snp_info{$index}->{"W"}->{"qual"}+$candidate_snp_info{$index}->{"C"}->{"qual"};
			if($allele){
				$genotype1 = $callRuleW{$allele}->{$call_w};
				$genotype2 = $callRuleC{$allele}->{$call_c};
				if($genotype1 eq "?" || $genotype2 eq "?"){
					$call = "?";
					$call = $genotype1 if($genotype2 eq "?");
					$call = $genotype2 if($genotype1 eq "?");
				}else{
					$call = $callRuleDs{$genotype1}->{$genotype2};
				}
			}else{
				$call = $callRuleDs{$call_w}->{$call_c};
				$call = 0 if($candidate_snp_info{$index}->{"W"}->{"depth"} < 5 || $candidate_snp_info{$index}->{"C"}->{"depth"} < 5);
			}
		}elsif($allele){
			if($call_w){
				$call = $callRuleW{$allele}->{$call_w};
				$call = "?" if($call eq "Y");
				$qual = $candidate_snp_info{$index}->{"W"}->{"qual"};
			}elsif($call_c){
				$call = $callRuleC{$allele}->{$call_c};
				$call = "?" if($call eq "R");
				$qual = $candidate_snp_info{$index}->{"C"}->{"qual"};
			}
		}
		$candidate_snp_info{$index}->{"call"}=$call if($call && $call ne "?");
		$candidate_snp_info{$index}->{"qual"}=$qual;
	}
}

sub load_dbsnp(){
	my $refSnpFile = shift;
	open(INFILE, "$refSnpFile") || die("Error opening $refSnpFile\n");
	while(my $line = <INFILE>){
		chomp($line);
		my @fields = split(/\t/, $line);
		my $index = $fields[1] . ":" . $fields[3];
		if($candidate_snp_info{$index}){
			next if(!$three2one{$fields[9]});
			next if($candidate_snp_info{$index}->{"rs"} && $candidate_snp_info{$index}->{"rs"} =~ m/rs/);
			$candidate_snp_info{$index}->{"rs"} = $fields[4];
			$candidate_snp_info{$index}->{"allele"} = $fields[6] eq "+" ? $three2one{$fields[9]} : $rcTable{$three2one{$fields[9]}};
		}
	}
	close(INFILE);
}

sub recalibrate{
	my ($str, $call, $allele, $count) = @_;
	my %rules;
	$allele = "N" if(!$allele);
	%rules = %hWRules if($str eq 'W');
	%rules = %hCRules if($str eq 'C'); 
	my @f = split ",", $count;
	my ($total, $sum_not, $second) = (0,0, "N");
	for(my $i = 0; $i < scalar(@f)-1; $i+=2){
		$total+=$f[$i+1];
		$sum_not+=$f[$i+1] if($f[$i] ne $call);
		$second = $f[$i] if($f[$i] ne $call);
	}
	my $geno = ord($call) < ord($second) ? $call . "/" . $second : $second . "/". $call;
	$geno = $three2one{$geno};
	if($total - $sum_not >= $total*0.9){
		my $score = est_binomial_phred($total, $total-$sum_not);
		if($sum_not == 0 || $score > 50){
			return $call;
		}elsif( scalar(@f) == 4 && $score <= 50 && $sum_not > 1){
			my $new_call = $rules{$allele}->{$call};
			if($allele eq "N"){
				$call = 0;
			}else{
				$geno = $callRuleW{$allele}->{$geno} if($str eq 'W');
				$geno = $callRuleC{$allele}->{$geno} if($str eq 'C');
				$call = $new_call if($new_call && $new_call eq $geno);
			}
			return $call;
		}else{
			return 0;
		}
	}elsif(scalar(@f) == 4){
		my $new_call = $rules{$allele}->{$call};
		if($allele eq "N"){
			$call = 0;
		}else{
			$geno = $callRuleW{$allele}->{$geno} if($str eq 'W');
			$geno = $callRuleC{$allele}->{$geno} if($str eq 'C');
			$call = $new_call if($new_call && $new_call eq $geno);
		}
		return $call;
	}else{
		return 0;
	}
}

sub est_binomial_phred{
	my ($n, $k) = @_;
	my $F = exp( (-16/$n)*($n/2-$k)**2 )/15;
	my $phred = 100000;
	$phred = -10*log10($F) if($F > 0);
	return $phred;
}
sub log10{
	my $n = shift;
	return log($n)/log(10);
}

main();
