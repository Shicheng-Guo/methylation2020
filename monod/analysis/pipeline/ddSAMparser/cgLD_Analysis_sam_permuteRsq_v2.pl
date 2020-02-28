#!/usr/bin/perl -w
# this script reads in a SAM format file with fillmd field, and output the rsq value for each pair of cg sites

use strict;
use warnings;
use Sort::Array qw/Sort_Table/;

my %cgTable;

my %RsqTable;
my %RsqPvalTable;

my %methylTable;
my $min_coverage = 10;
my $qual_base = 64;
my $min_rsq = 0.3;

my %alignment;
my $deletion_string = "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD";

my %rcTable;
$rcTable{'A'}='T';
$rcTable{'T'}='A';
$rcTable{'G'}='C';
$rcTable{'C'}='G';
$rcTable{'N'}='N';
$rcTable{'D'}='D';

my @palette=("0,240,0", "30,210,0","60,180,0","90,150,0","120,120,0","150,90,0","180,60,0","210,0,0");

sub main{
	printUsage() if(!$ARGV[0] || !$ARGV[1]);
	readCGTable();
	while(my $line = <STDIN>){
		chomp($line);
		next if($line !~ m/MD:Z/);
		my %alignment = getAlignmentInfo($line);
		getReference(\%alignment);
		getCGPos(\%alignment);
		next if(!$alignment{"cgPos"});
		my @cgPos = @{$alignment{"cgPos"}};
		my %indices;
		foreach my $pos_call (@cgPos){
			my ($pos, $call) = split ":", $pos_call;
			if($indices{$pos}){
				# pos has already been assigned?
				die("Getting ambiguous call for position $pos for:\n$line\n");
				delete $indices{$pos};
			}else{
				#$methylTable{$alignment{"chr"}.":".$pos}->{$call}++;
				next if($call !~ /[CT]/);
				$indices{$pos} = $call;
			}
		}
		my @ordered_pos = sort {$a<=>$b} keys %indices;
		next if(scalar(@ordered_pos) < 2);
		#print join(",", @ordered_pos), "\n";
		for(my $i = 0; $i < scalar(@ordered_pos); $i++){
			for(my $j = $i+1; $j < scalar(@ordered_pos); $j++){
				my $value_key = $alignment{"chr"}.":".$ordered_pos[$i]."-".$ordered_pos[$j];
				my $value_call = $indices{$ordered_pos[$i]}.$indices{$ordered_pos[$j]};
				#$value_call =~ tr/T/U/;
				#$value_call =~ tr/C/M/;
				my ($mm, $mu, $um, $uu) = (0,0,0,0);
				my ($rsq, $pval, $counts, $totalHaps) = ("NA", "NA", "NA", 0);
				if($RsqTable{$value_key}){
					($rsq, $pval, $counts, $totalHaps) = split "\t", $RsqTable{$value_key};
					($mm, $mu, $um, $uu) = split ":", $counts;
				}
				$mm++ if($value_call eq "CC");
				$mu++ if($value_call eq "CT");
				$um++ if($value_call eq "TC");
				$uu++ if($value_call eq "TT");
				$totalHaps++;
				$counts = "$mm:$mu:$um:$uu";
				if($totalHaps >= $min_coverage){
					# uniform site is where one cpg have maf of 0
					my $mafA = ($um + $uu) < ($mm + $mu) ? $um + $uu : $mm + $mu;
					my $mafB = ($mu + $uu) < ($mm + $um) ? $mu + $uu : $mm + $um;
					if($mafA ne 0 and $mafB ne 0){
						($rsq, $pval) = updateRsq($counts);
						$RsqTable{$value_key} = "$rsq\t$pval\t$counts\t$totalHaps";
					}else{
						$RsqTable{$value_key} = "0\tMAF_0\t$counts\t$totalHaps";
					}
				}else{
					$RsqTable{$value_key} = "$rsq\t$pval\t$counts\t$totalHaps";
				}
			}
		}
		undef(@ordered_pos);
		undef(%indices);
		undef(%alignment);
	}
	outputRsqTable();
	undef(%RsqTable);
	undef(%RsqPvalTable);
	undef(%cgTable);
	#printBEDTable();
	undef(%methylTable);
}

sub printBEDTable{
	my $sampleID = $ARGV[0];
	my $file_name = $sampleID . ".BED.txt";
	open(OUT, ">$file_name") || die("Error writing $file_name\n");
	print OUT "track name=\"", $sampleID, "\" description=\"$sampleID Methylation level\" visibility=2 useScore=1 itemRgb=\"On\"\n"; 
	foreach my $chr_pos (keys %methylTable){
		my ($chr, $pos) = split /:/, $chr_pos;
		my %allele_counts;
		my $CT_count = 0;
		my $total_depth = 0;
		foreach my $call (keys %{$methylTable{$chr_pos}}){
			my $depth = $methylTable{$chr_pos}->{$call};
			$allele_counts{$call} = $depth;
			$CT_count += $depth if($call =~ /[CT]/);
			$total_depth += $depth;
		}
		next if($total_depth < $min_coverage);
		next if(!$CT_count || $CT_count/$total_depth < 0.9);
		my $m = $allele_counts{"C"} ? $allele_counts{"C"} : 0;
		my $methylLevel = $m/$CT_count;
		print OUT "$chr\t", $pos-1, "\t$pos\t\'$m/$CT_count\'\t$CT_count\t+\t", $pos-1, "\t$pos\t", $palette[int($methylLevel*8-0.0001)], "\n";
	}
	close(OUT);
}

sub readCGTable{
	open(IN, "$ARGV[1]") || die("Error opening $ARGV[1]\n");
	while(my $line = <IN>){
		chomp($line);
		my @tmp = split /\t/, $line;
		$tmp[0] =~ s/:W//;
		$cgTable{$tmp[0].":".$tmp[1]} = 1;		
	}
	close(IN);
}

sub getReference{
	my $alignment = shift;
	my $cigar = $alignment->{"cigar"};
	my $query = $alignment->{"seq"};
	my $qual = $alignment->{"qual"};
	my $ref_seq;
	my $ref_qual;
	# Adapted from Ben Langmead (BtlBio::Alignment:Util.pm)
	# CIGAR fields are always in pairs
	my $i = 0;
	my $j = 0;
	my $nm_i = 0;
	my $nm_d = 0;
	while($i < length($cigar)){
		substr($cigar, $i) =~ /^([0-9]+)/;
		defined($1) || die("Could not parse number at pos $i: '$cigar'");
		my $runlen = $1;
		$i += length($1);
		$i < length($cigar) || die("Bad cigar string : '$cigar'");
		my $op = substr($cigar, $i, 1);
		defined($op) || die("count not parse operation at pos $i: '$cigar'");
		die("Could not understand $op: '$cigar'") if($op !~ m/[MX=DIS]/);
		$op =~ s/[X=]/M/g;
		my ($clip_s, $clip_q);
		if($op eq "M" || $op eq "I" || $op eq "S"){
			$clip_s = substr($query, $j, $runlen);
			$clip_q = substr($qual, $j, $runlen);
			$clip_s =~ s/[ATGCatgc]/I/g if($op eq "I");
			$nm_i += $runlen if($op eq "I");
			$j += $runlen;
		}else{
			#deletion from reference
			$nm_d += $runlen;
			length($deletion_string) > $runlen || die("deletion is too long at $runlen: '$cigar'");
			$clip_s = substr($deletion_string, 0, $runlen);
			$clip_q = substr($deletion_string, 0, $runlen);
		}
		$i++;
		$ref_seq = $ref_seq . $clip_s if($op =~ m/[MD]/);
		$ref_qual = $ref_qual . $clip_q if($op =~ m/[MD]/);
	}
	#print "Ref_match_seq: $ref_seq\n";
	#print "Ref_match_qual: $ref_qual\n";

	$alignment->{"ref_match_seq"} = $ref_seq;
	$alignment->{"ref_match_qual"} = $ref_qual;
	$alignment->{"ref_nm_i"} = $nm_i;
	$alignment->{"ref_nm_d"} = $nm_d;
	$alignment->{"tlen"} = length($ref_seq) - $nm_i if($ref_seq);
}

sub getCGPos{
	my $alignment = shift;
	my @mismatches;
	my $mdz = $alignment->{"mdz"};
	my $ref_seq = uc($alignment->{"ref_match_seq"});
	my $ref_qual = $alignment->{"ref_match_qual"};
	return @mismatches if(!$alignment->{"ref_match_seq"});
	return @mismatches if(!$alignment->{"ref_match_qual"});
	#print $ref_seq, "\n", $ref_qual, "\n";
	my $chr = $alignment->{"chr"};
	my $position = $alignment->{"pos"};
	my @tmp = split /(\d+)/, $mdz;
	my @bases = split "", $ref_seq;
	my @quals = split "", $ref_qual;
	my $j = 0;
	my $relative_pos = 0;
	for(my $i = 1; $i< scalar(@tmp); $i++){
		my $op = $tmp[$i];
		if($op =~ s/^\^//){
			$j+=length($op);
			$relative_pos+=length($op);
			#print $op, "\n";
		}elsif($op =~ /^[+-]?\d+$/){
			for(my $k = 0; $k < $op; $k++){
				# obtain calls at the current position in ref_seq
				my $b = $bases[$j];
				my $q = $quals[$j];
				my $actual_pos = $position + $relative_pos;
				$j++; # j position moves along ref_seq
				if($b eq "I"){ $k--; next;} # if the position is an insertion to the reference, we need to go back 1
				$relative_pos++; # relative_position on get incremented when the position is not an insertion to the reference
				next if(ord($q) - $qual_base < 5); # do not make a call if the quality is bad
				#print "SM: $relative_pos\t$op\t$b\t$q\n";
				if($alignment->{"orientation"} eq "W"){
					push(@{$alignment->{"cgPos"}}, "$actual_pos:$b") if($cgTable{"$chr:$actual_pos"});
				}elsif($alignment->{"orientation"} eq "C"){
					$actual_pos--;
					die("Unrecognizable base $b\n") if(!$rcTable{$b});
					my $rc_b = $rcTable{$b};
					push(@{$alignment->{"cgPos"}}, "$actual_pos:$rc_b") if($cgTable{"$chr:$actual_pos"});
				}
			}
		}elsif($op =~ /[ATGCatgc]/){
			# obtain calls at the current position in ref_seq
			my $b = $bases[$j];
			my $q = $quals[$j];
			my $actual_pos = $position + $relative_pos;
			# make the necessary incrementations in j and relative positions
			$relative_pos++;
			$j++;
			next if(ord($q) - $qual_base < 5); # do not make a call if the quality is bad
			#print "NM: $relative_pos\t$op\t$b\t$q\n";
			if($alignment->{"orientation"} eq "W"){
				push(@{$alignment->{"cgPos"}}, "$actual_pos:$b") if($cgTable{"$chr:$actual_pos"});
			}elsif($alignment->{"orientation"} eq "C"){
				$actual_pos--;	
				die("Unrecognizable base $b\n") if(!$rcTable{$b});
				my $rc_b = $rcTable{$b};
				push(@{$alignment->{"cgPos"}}, "$actual_pos:$rc_b") if($cgTable{"$chr:$actual_pos"});
			}
		}
	}
	undef @bases;
	undef @quals;
	undef @tmp;
	return @mismatches;
}

sub getAlignmentInfo{
	my $line = shift;
	my @tmp = split /\t/, $line;
	my %alignment;
	$alignment{"name"} = $tmp[0];
	$alignment{"orientation"} = $tmp[1] & 16 ? "C":"W";
	$alignment{"chr"} = $tmp[2];
	$alignment{"pos"} = $tmp[3];
	$alignment{"cigar"} = $tmp[5];
	$alignment{"seq"} = $tmp[9];
	$alignment{"qual"} = $tmp[10];
	$alignment{"ndz"} = $tmp[11];
	$alignment{"mdz"} = $tmp[12];
	$alignment{"ref_match_seq"} = "NA";
	$alignment{"ref_match_qual"} = "NA";
	$alignment{"ref_nm_i"} = 0;
	$alignment{"ref_nm_d"} = 0;
	$alignment{"tlen"} = 0;
	$alignment{"cgPos"} = ();
	return %alignment;
}

sub outputRsqTable{
	my $chr;
	my $dist;
	open(RSQ_OUT, ">$ARGV[0]") || die("error writing to file $ARGV[0]\n");
	print RSQ_OUT "index\tRsq\tmm:mu:um:uu\ttotalHaps\n";
	foreach my $pairs(keys %RsqTable){
		my ($rsq, $pval, $counts, $totalHaps) = split "\t", $RsqTable{$pairs};
		next if($totalHaps < $min_coverage);
		print RSQ_OUT $pairs, "\t", $RsqTable{$pairs},"\n";	
	}
	close(RSQ_OUT);
}

sub updateRsq{
	my $counts = shift;
	my ($mm, $mu, $um, $uu) = split ":", $counts;
	my $total = $mm+$mu+$um+$uu;
	#calculate p1 = prob(1M), p2=prob(1U), q1=prob(2M), q2=prob(2U)
	#calculate D = freq(1M/2M) - prob(1M)*prob(2M)
	#calculate r^2 = D^2/(p1*p2*q1*q2)
	my ($p1, $q1, $p2, $q2) = ($mm+$mu, $mm+$um, $um+$uu, $mu+$uu);
	$p1 /=$total;
	$q1 /=$total;
	$p2 /=$total;
	$q2 /=$total;
	$mm /=$total;
	$mu /=$total;
	$um /=$total;
	$uu /=$total;
	my $D = $mm*$uu - $um*$mu;
	my $rsq = 0;
	$rsq = ($D*$D)/($p1*$p2*$q1*$q2);
	my $out_rsq = sprintf("%4.3f", $rsq);
	my $out_pval;
	my ($mm_s, $mu_s, $um_s, $uu_s) = (int($min_coverage*$mm), int($min_coverage*$mu), int($min_coverage*$um), int($min_coverage*$uu));
	my $counts_s = "$mm_s:$mu_s:$um_s:$uu_s";
	if(defined($RsqPvalTable{$counts_s})){
		return $out_rsq, $RsqPvalTable{$counts_s};
	}
	# begins permutation:
	my $num_perm = 100;
	$p1 = sprintf("%.3f", $p1);
	$q1 = sprintf("%.3f", $q1);
	my $ptotal = $min_coverage;
	my $greater_eq_than = 0;
	for(my $i = 0; $i < $num_perm; $i++){
		my ($r_mm, $r_mu, $r_um, $r_uu) = (0,0,0,0);
		for(my $j = 0; $j < $ptotal; $j++){
			my $a = int(rand(1001))/1000;
			my $b = int(rand(1001))/1000;
			$r_mm++ if($a <= $p1 and $b <= $q1);
			$r_mu++ if($a <= $p1 and $b > $q1);
			$r_um++ if($a > $p1 and $b <= $q1);
			$r_uu++ if($a > $p1 and $b > $q1);
		}
		$r_mm /= $ptotal;
		$r_mu /= $ptotal;
		$r_um /= $ptotal;
		$r_uu /= $ptotal;
		my $r_D = $r_mm*$r_uu - $r_um*$r_mu;
		my $r_rsq = ($r_D*$r_D)/($p1*$p2*$q1*$q2);
		$greater_eq_than++ if($r_rsq >= $rsq);
	}
	$out_pval = sprintf("%4.3f", $greater_eq_than/$num_perm);
	$RsqPvalTable{$counts_s} = $out_pval;
	return $out_rsq, $out_pval;
}

sub printUsage{
	print " Usage: \n";
	print " ./cgLD_Analysis_sam_permuteRsq_v2.pl [rsq_out_put_name] [cpg position list] < fillmd.sam\n";
	exit 0;
}

main();
