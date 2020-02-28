# This software is Copyright © 2017 The Regents of the University of California. All Rights Reserved.
#  
# Permission to copy, modify, and distribute this software and its documentation for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
#  
# Permission to make commercial use of this software may be obtained by contacting:
# Office of Innovation and Commercialization
# 9500 Gilman Drive, Mail Code 0910
# University of California
# La Jolla, CA 92093-0910
# (858) 534-5815
# kzhang@eng.ucsd.edu
 
# This software program and documentation are copyrighted by The Regents of the University of California. The software program and documentation are supplied "as is", without any accompanying services from The Regents. The Regents does not warrant that the operation of the program will be uninterrupted or error-free. The end-user understands that the program was developed for research purposes and is advised not to rely exclusively on the program for any reason.
#  
# IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
# CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
# EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF
# CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
# MODIFICATIONS.


#!/home/shg047/perl5/perlbrew/perls/perl-5.22.0/bin/perl -w
# this script reads in a haploInfo file and output the rsq value for each pair of cg sites
# Transfer Bam to Fastq with samtools command
# Run the script to the Bam directory
# Contact: Shicheng Guo <shihcheng.guo@gmail.com>
# Version 1.3
# Update: 2016-02-25

use strict;
use strict;
use warnings;
use Sort::Array qw/Sort_Table/;

my %cgTable;

my %RsqTable;
my %RsqPvalTable;
my %posList;

my %methylTable;
my $min_coverage = 2;
my $min_rsq = 0.1;
my $pvalue_cutoff = 0.1;
my $distance_threshold = 200;

my %rcTable;
$rcTable{'A'}='T';
$rcTable{'T'}='A';
$rcTable{'G'}='C';
$rcTable{'C'}='G';
$rcTable{'N'}='N';
$rcTable{'D'}='D';

printUsage() if(!$ARGV[0]);
my ($target_chr, $target_start, $target_end) = split /[:-]/, $ARGV[1];

sub main{
	while(my $line = <STDIN>){
		next if $line=~/^\s+$/;
         	chomp($line);
		my @tmp =  split /\s+/, $line;
                next if scalar(@tmp)<4;
                my ($chr, $start, $end) = split /[:-]/, $tmp[0];
		next if($chr ne $target_chr || $start > $target_end || $end < $target_start);
		my ($hapString, $hapCount) = ($tmp[1], $tmp[2]);
		my @cgPos = split ",", $tmp[3];
		next if(length($tmp[1]) ne scalar(@cgPos));
		my %indices;
		for(my $i = 0; $i < scalar(@cgPos); $i++){
			my $pos = $cgPos[$i];
			my $call = substr($hapString,$i,1);
			if($indices{$pos}){
				# pos has already been assigned?
				die("Getting ambiguous call for position $pos for:\n$line\n");
				delete $indices{$pos};
			}else{
				next if($call !~ /[CT]/);
				$indices{$pos} = $call;
			}
		}
		my @ordered_pos = sort {$a<=>$b} keys %indices;
		next if(scalar(@ordered_pos) < 2);
		#print join(",", @ordered_pos), "\n";
		for(my $i = 0; $i < scalar(@ordered_pos); $i++){
			for(my $j = $i+1; $j < scalar(@ordered_pos); $j++){
				my $value_key = $chr.":".$ordered_pos[$i]."-".$ordered_pos[$j];
				my $value_call = $indices{$ordered_pos[$i]}.$indices{$ordered_pos[$j]};
				#$value_call =~ tr/T/U/;
				#$value_call =~ tr/C/M/;
				my ($mm, $mu, $um, $uu) = (0,0,0,0);
				my ($rsq, $pval, $counts, $totalHaps) = ("NA", 0, "NA", 0);
				if($RsqTable{$value_key}){
					($rsq, $pval, $counts,$totalHaps) = split "\t", $RsqTable{$value_key};
					($mm, $mu, $um, $uu) = split ":", $counts;
				}
				$mm+=$hapCount if($value_call eq "CC");
				$mu+=$hapCount if($value_call eq "CT");
				$um+=$hapCount if($value_call eq "TC");
				$uu+=$hapCount if($value_call eq "TT");
				$totalHaps+=$hapCount;
				$counts = "$mm:$mu:$um:$uu";
				if($totalHaps >= $min_coverage){
					# uniform site is where one cpg have maf of 0
					my $mafA = ($um + $uu) < ($mu + $mm) ? $um + $uu : $mu + $mm;
					my $mafB = ($mu + $uu) < ($um + $mm) ? $mu + $uu : $um + $mm;
					if($mafA ne 0 and $mafB ne 0){
						($rsq, $pval) = updateRsq($counts);
						$RsqTable{$value_key} = "$rsq\t$pval\t$counts\t$totalHaps";
					}else{
						$rsq = 1 if($mafA + $mafB == 0);
						$rsq = ($mu+$um)/$mm if($uu == 0 && $mafA + $mafB > 0);
						$rsq = ($mu+$um)/$uu if($mm == 0 && $mafA + $mafB > 0);
						if($rsq ne "NA"){
							$rsq = 1/$rsq if($rsq > 1);
							$rsq = $rsq > (1-$rsq) ? $rsq : 1 - $rsq;
							my $out_rsq = sprintf("%4.3f", $rsq);
							$RsqTable{$value_key} = "$out_rsq\tMAF_0\t$counts\t$totalHaps";
						}else{
							$RsqTable{$value_key} = "NA\tMAF_0\t$counts\t$totalHaps";
						}
					}
				}else{
					$RsqTable{$value_key} = "$rsq\t$pval\t$counts\t$totalHaps";
				}		
			}
		}
		undef(@ordered_pos);
		undef(%indices);
	}
	outputRsqTable();
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

sub reportBlocks{
	my $filename = shift;
	my $h_blocks = shift;
	my %blocktable = %{$h_blocks};
	open(OUTFILE, ">$filename") || die("Error in opening file $filename\n");
	print OUTFILE "block id\tstart locus\tend locus\n";
	foreach my $id (sort keys(%blocktable)) {
		next if($blocktable{$id}->{'end'} - $blocktable{$id}->{'start'} < 100);
		print OUTFILE $id, "\t", $blocktable{$id}->{'start'}, "\t", $blocktable{$id}->{'end'}, "\n";
	}
	close(OUTFILE);
}

sub outputRsqTable{
	# open(RSQ_OUT, ">$ARGV[0]") || die("error writing to file $ARGV[0]\n");
	my %table;
	#print RSQ_OUT "index\tRsq\tpval\tmm:mu:um:uu\ttotalHaps\n";
	foreach my $pairs(keys %RsqTable){
		my ($chr, $start, $end) = split /[:-]/, $pairs;
		my ($rsq, $pval, $counts, $totalHaps) = split "\t", $RsqTable{$pairs};
		next if($totalHaps < $min_coverage);
		$posList{$chr}->{$start} = 1;
		$posList{$chr}->{$end} = 1;
		#print RSQ_OUT $pairs, "\t", $RsqTable{$pairs},"\n";
		$table{$start}->{$end} = $rsq;
		$table{$end}->{$start} = $rsq;
	}
	foreach my $chr (keys %posList){
                open(RSQ_OUT, ">$ARGV[0].$chr.rsq") || die("error writing to file $ARGV[0].chr.rsq\n");
		my @sorted_list = sort {$a<=>$b} keys %{$posList{$chr}};
		print RSQ_OUT "\t", join("\t", @sorted_list), "\n"; 
		for(my $i = 0; $i < scalar(@sorted_list); $i++){
			print RSQ_OUT $sorted_list[$i];
			for(my $j = 0; $j < scalar(@sorted_list); $j++){
				if(!$table{$sorted_list[$i]}->{$sorted_list[$j]}){
					print RSQ_OUT "\tNA";
					next;
				}
				print RSQ_OUT "\t", $table{$sorted_list[$i]}->{$sorted_list[$j]};
			}
			print RSQ_OUT "\n";
		}
		close(RSQ_OUT);
	}
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
	my $abs_Dprime;
	my $Dmax = $D > 0 ? ($p1*$q1 < $p2*$q2 ? $p1*$q1:$p2*$q2) : ($p1*$q2 < $p2*$q1 ? $p1*$q2:$p2*$q1);
	my $rsq = 0;
	if($D == 0){
		$abs_Dprime = ($Dmax == 0.0) ? 1.0 : 0.0;
		$rsq = ($p1*$p2*$q1*$q2) == 0 ? 1.0 : 0.0;
	}else{
		$abs_Dprime = abs($D/$Dmax);
		$rsq = ($D*$D)/($p1*$p2*$q1*$q2);
	}
	my $out_rsq = sprintf("%4.3f", $rsq);
	my $out_pval;
	my ($mm_s, $mu_s, $um_s, $uu_s) = (int($min_coverage*$mm), int($min_coverage*$mu), int($min_coverage*$um), int($min_coverage*$uu));
	my $counts_s = "$mm_s:$mu_s:$um_s:$uu_s";
	if(defined($RsqPvalTable{$counts_s})){
		return $out_rsq, $RsqPvalTable{$counts_s};
	}
	# begins permutation:
	my $num_perm = 500;
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

sub lookupLD{
	my $chr = shift;
	my $pos1 = shift; # first cpg position
	my $pos2 = shift; # last cpg position
	my $rsq;
	if(!$RsqTable{$chr.":".$pos1."-".$pos2}){
		$rsq = "NA"; # no evidence for LD so split block
	}else{
		my ($rsq_out, $pval, $allele_count, $totalHaps) = split "\t", $RsqTable{$chr.":".$pos1."-".$pos2};
		$rsq = $rsq_out;
		$rsq = 0 if($pval ne "MAF_0" and $pval > $pvalue_cutoff);
	}
	return $rsq;
}

sub getBlocksGreedy{
	my $sorted_list = shift;
	my $chr = shift;
	my $start = shift;
	my $end = shift;
	my $threshold = shift;
	my $i = $start;
	my @blocks;
	my %blockSetTable;
	print "Total nummber of CpGs is $end for $chr\n";
	while($i <= $end){
		# continue adding CpGs to set if there is another CpG within distance threshold.
		while(($i < $end) and @{$sorted_list}[$i+1] - @{$sorted_list}[$i] < $distance_threshold){
			$i++;
		}
		print $i, "\n";
		push(@blocks, getAllBlockInRegion($sorted_list, $chr, $start, $i, $threshold));
		$i++;
		$start = $i;
	}
	my @sortedBlocks = Sort_Table(
			cols		=> '2',
		        field		=> '1',
			sorting  	=> 'ascending',
			structure	=> 'csv',
			separator	=> '\:',
			data		=> \@blocks,
	);
	for(my $i=0; $i<scalar(@blocks); $i++){
		my $id = sprintf("B%05d", $i+1);
		my @words = split(/:/, $sortedBlocks[$i]);
		$blockSetTable{$id}->{'start'} = $words[0];
		$blockSetTable{$id}->{'end'} = $words[1];
	}
	return %blockSetTable;
}

sub getAllBlockInRegion{
	my $sorted_list = shift;
	my $chr = shift;
	my $start = shift;
	my $end = shift;
	my $threshold = shift;
	my @blocks;
	if($start == $end){
		push(@blocks, @{$sorted_list}[$start] . ":" . @{$sorted_list}[$end]);
		return @blocks;
	}
	my ($block_start, $block_end) = findMaxBlockInRegion($sorted_list, $chr, $start, $end, $threshold);
	push(@blocks, @{$sorted_list}[$block_start] . ":" . @{$sorted_list}[$block_end]) if($block_end - $block_start > 2);
	# if the largest block is not all inclusive:
	if($block_start > $start){
		my @sub_blocks = getAllBlockInRegion($sorted_list, $chr, $start, $block_start-1, $threshold);
		push(@blocks, @sub_blocks);
	}
	if($block_end < $end){
		my @sub_blocks = getAllBlockInRegion($sorted_list, $chr, $block_end+1, $end, $threshold);
		push(@blocks, @sub_blocks);
	}
	return @blocks;
}

sub findMaxBlockInRegion(){
	my $sorted_list = shift;
	my $chr = shift;
	my $start = shift;
	my $end = shift;
	my $threshold = shift;
	my $max_block_start=$start;
	my $max_block_end=$start;
	return ($start, $end) if($start == $end);
	for(my $size = $end-$start+1; $size >1; $size--){ # try to get the largest size first
                my $good_block=0;
                for(my $i= $start; $i<= $end-$size+1; $i++){ # allow every pos to be starting point
                        $good_block = 1;
                        for(my $j= $i; $j<$i+$size; $j++){ # every pos to every other pos within block must be in LD!
                                for(my $k = $j+1; $k<$i+$size; $k++){
                                        if(lookupLD($chr, @{$sorted_list}[$j], @{$sorted_list}[$k]) ne "NA") {
                                                $good_block = 0 if(lookupLD($chr, @{$sorted_list}[$j], @{$sorted_list}[$k]) < $threshold);
                                                last;
                                        }
                                }
                        }
                        if($good_block){
                                $max_block_start = $i;
                                $max_block_end = $i + $size -1;
                                last;
                        }
                }
                last if($good_block);
        }
	return ($max_block_start,  $max_block_end);
}

sub printUsage{
        print " Usage:\n";
		print " perl ~/monod/bin/$0 <rlt_prefix> <chr:start-end> < haploinfo.input.txt\n";
		print " For example:\n perl $0 rlt chr10:10000873-10001472 < haploinfo.input.txt\n\n";
		print "-----------------------------------------------------------------------------\n";
		print " haploinfo.input.txt format:\n";
		print " chr10:10000873-10001472        CCC     1       10001056,10001082,10001168\n\n";
    	print "-----------------------------------------------------------------------------\n";
        exit 0;
}


main();

