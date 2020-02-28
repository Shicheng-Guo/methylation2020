#!/usr/bin/perl -w
use strict;

# CONSTANTS
my $smoothing = 0.001;
my $minHapLen = 4;
my $minimal_mean_mC_diff = 0.2;

# optimizable parameters:
my $min_hap_counts = $ARGV[3] ? $ARGV[3] : 10;
my $mean_error_rate = $ARGV[4] ? $ARGV[4] : 0.03;
my $min_evidence = $ARGV[5] ? $ARGV[5] : 0.01;

# Data structures:
my %read_counts_by_chr;
my %cgTable;
my %binompTable;
# %hapTable => haps => hapCount
#           => totalHaps
#	    => cpgPositions

# USER INPUTS
my $ref_hap_file_A = $ARGV[0];
my $ref_hap_file_B = $ARGV[1];
my $mix_hap_file = $ARGV[2];

#$ref_hap_file_A =~ s/.haploInfo.txt//g;
#$ref_hap_file_B =~ s/.haploInfo.txt//g;
#$mix_hap_file =~ s/.haploInfo.txt//g;

sub main{
	
	my %hapTableA = load_data_file($ref_hap_file_A);
	my %hapTableB = load_data_file($ref_hap_file_B);
	my %hapTableMix = load_data_file($mix_hap_file);

	my ($probe_used, $hap_used, $total_assigned_to_A, $total_assigned_to_B, $total_unassigned) = (0,0,0,0,0);
	my ($total_log_score_to_A, $total_log_score_to_B) = (0,0);

	open(OUT_DATA, ">MinorHaps.txt") || die("Error writing to MinorHaps.txt\n");
	foreach my $probeID (keys %hapTableMix){
		my %matchedHapTable;
		$matchedHapTable{"totalHap"}->{"A"} = 0;
		$matchedHapTable{"totalHap"}->{"B"} = 0;
		my ($chr,$chr_start,$chr_end) = split /[:-]/, $probeID;
		my @cpgPositions = @{$hapTableMix{$probeID}->{"cpgPositions"}};
		next if(scalar(@cpgPositions) < $minHapLen || $chr =~ /[XYM]/);
		my %cgOverlaps;
		#print join(":", @cpgPositions), "\n";
		foreach my $pos (@cpgPositions){
			next if(!$cgTable{$chr.":".$pos});
			my @candidates = keys %{$cgTable{$chr.":".$pos}};
			foreach my $canProbeID (@candidates){
				#print $canProbeID, "\n";
				$cgOverlaps{$canProbeID}->{"A"}++ if($hapTableA{$canProbeID});
				$cgOverlaps{$canProbeID}->{"B"}++ if($hapTableB{$canProbeID});
			}
		}
		#print $probeID,"\n";
		foreach my $canProbeID (keys %cgOverlaps){
			if($cgOverlaps{$canProbeID}->{"A"} && $cgOverlaps{$canProbeID}->{"A"} >= $minHapLen){
				#print $canProbeID, " in A\n";
				my @cpgPositions_in_A = @{$hapTableA{$canProbeID}->{"cpgPositions"}};
				foreach my $hapString (keys %{$hapTableA{$canProbeID}->{"haps"}}){
					# trim/pad so that hapString aligns with the mix hapString
					my $hapString_substr = align_hapString_to_pos(\@cpgPositions, \@cpgPositions_in_A, $hapString);
					#print "\t", $hapString_substr, " in A \n";
					$matchedHapTable{"hapCounts"}->{"A"}->{$hapString_substr} += $hapTableA{$canProbeID}->{"haps"}->{$hapString};
					$matchedHapTable{"totalHap"}->{"A"} += $hapTableA{$canProbeID}->{"haps"}->{$hapString};
				}
			}
			if($cgOverlaps{$canProbeID}->{"B"} && $cgOverlaps{$canProbeID}->{"B"} >= $minHapLen){
				#print $canProbeID, " in B\n";
				my @cpgPositions_in_B = @{$hapTableB{$canProbeID}->{"cpgPositions"}};
				foreach my $hapString (keys %{$hapTableB{$canProbeID}->{"haps"}}){
					# trim/pad
					my $hapString_substr = align_hapString_to_pos(\@cpgPositions, \@cpgPositions_in_B, $hapString);
					#print "\t", $hapString_substr, " in B \n";
					$matchedHapTable{"hapCounts"}->{"B"}->{$hapString_substr} += $hapTableB{$canProbeID}->{"haps"}->{$hapString};
					$matchedHapTable{"totalHap"}->{"B"} += $hapTableB{$canProbeID}->{"haps"}->{$hapString};
				}
			}
		}
		next if( $matchedHapTable{"totalHap"}->{"A"} < $min_hap_counts || $matchedHapTable{"totalHap"}->{"B"} < $min_hap_counts);
		my $refA_mean_methylation = sprintf("%4.3f", calc_mean_methylation(\%{$matchedHapTable{"hapCounts"}->{"A"}}));
		my $refB_mean_methylation = sprintf("%4.3f", calc_mean_methylation(\%{$matchedHapTable{"hapCounts"}->{"B"}}));
		next if(abs($refB_mean_methylation-$refA_mean_methylation)<$minimal_mean_mC_diff);
		$probe_used++;
		$hap_used += $hapTableMix{$probeID}->{"totalHap"};
		$read_counts_by_chr{$chr}->{"refA"} = 0 if(!$read_counts_by_chr{$chr}->{"refA"});
		$read_counts_by_chr{$chr}->{"refB"} = 0 if(!$read_counts_by_chr{$chr}->{"refB"});

		### ASSIGN HAPLOTYPES ###
		my ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned, $score_A, $score_B) = 
			assign_haps_by_similarity(\%{$matchedHapTable{"hapCounts"}->{"A"}}, \%{$matchedHapTable{"hapCounts"}->{"B"}}, \%{$hapTableMix{$probeID}->{"haps"}});
		#	assign_haps_by_similarity_NB(\%{$matchedHapTable{"hapCounts"}->{"A"}}, \%{$matchedHapTable{"hapCounts"}->{"B"}}, \%{$hapTableMix{$probeID}->{"haps"}});

		print OUT_DATA "$probeID\t$haps_assigned_to_B\t$haps_unassigned\t", $hapTableMix{$probeID}->{"totalHap"}, "\n";
		$read_counts_by_chr{$chr}->{"refA"}+=$haps_assigned_to_A;
		$read_counts_by_chr{$chr}->{"refB"}+=$haps_assigned_to_B;
		$read_counts_by_chr{$chr}->{"unassigned"}+=$haps_unassigned;
		$total_log_score_to_A += $score_A;
		$total_log_score_to_B += $score_B;
		undef(%matchedHapTable);
	}
	print "Mix=$mix_hap_file Probe_used=$probe_used Hap_used=$hap_used ";
	#print "CHROM\tRefA\tRefB\tUnassigned\n";
	foreach my $chr (keys(%read_counts_by_chr)){
		$read_counts_by_chr{$chr}->{"refA"}=0 if(!$read_counts_by_chr{$chr}->{"refA"});
		$read_counts_by_chr{$chr}->{"refB"}=0 if(!$read_counts_by_chr{$chr}->{"refB"});
		$read_counts_by_chr{$chr}->{"unassigned"}=0 if(!$read_counts_by_chr{$chr}->{"unassigned"});
		#print $chr, "\t", int($read_counts_by_chr{$chr}->{"refA"}), 
		#	"\t", int($read_counts_by_chr{$chr}->{"refB"}),
		#	"\t", int($read_counts_by_chr{$chr}->{"unassigned"}), "\n";
		$total_assigned_to_A+=$read_counts_by_chr{$chr}->{"refA"};
		$total_assigned_to_B+=$read_counts_by_chr{$chr}->{"refB"};
		$total_unassigned+=$read_counts_by_chr{$chr}->{"unassigned"}
	}
	my ($percent_A, $percent_B) = (0,0);
	if($total_assigned_to_A + $total_assigned_to_B > 0){
		$percent_A = sprintf("%4.3f", 100*$total_assigned_to_A/($total_assigned_to_A+$total_assigned_to_B));
		$percent_B = sprintf("%4.3f", 100*$total_assigned_to_B/($total_assigned_to_A+$total_assigned_to_B));
	}
	$total_assigned_to_A = sprintf("%.3f", $total_assigned_to_A);
	$total_assigned_to_B = sprintf("%.3f", $total_assigned_to_B);
	$total_log_score_to_A = sprintf("%.3f", $total_log_score_to_A);
	$total_log_score_to_B = sprintf("%.3f", $total_log_score_to_B);
	print "RefA=$total_assigned_to_A RefB=$total_assigned_to_B Unassigned=$total_unassigned ";
	print "$ref_hap_file_A=$percent_A (%), $total_log_score_to_A $ref_hap_file_B=$percent_B (%), $total_log_score_to_B\n";
	close(OUT_DATA);

	undef(%hapTableA);
	undef(%hapTableB);
	undef(%hapTableMix);
	undef(%read_counts_by_chr);
}

#align hapstring to a reference position list
sub align_hapString_to_pos{
	my $pos_list_ref = shift;
	my $pos_list_cur = shift;
	my $hapString = shift;
	my @bases = split "", $hapString;
	my %cgPosList;
	$cgPosList{@$pos_list_ref[$_]}->{"ref"} = $_ for (0..scalar(@$pos_list_ref)-1);
	$cgPosList{@$pos_list_cur[$_]}->{"cur"} = $_ for (0..scalar(@$pos_list_cur)-1);
	my $hapString_substr;
	foreach my $pos (sort {$a<=>$b} keys %cgPosList){
		if(!defined($cgPosList{$pos}->{"ref"})){ #not in ref
			# dont do anything
		}elsif(defined($cgPosList{$pos}->{"cur"})){ #in ref and in cur
			$hapString_substr = $hapString_substr . $bases[$cgPosList{$pos}->{"cur"}];
		}else{ # in ref and not in cur, pad
			$hapString_substr = $hapString_substr . "N";
		}
	}
	return $hapString_substr;
}

#assign haplotype by NB
sub assign_haps_by_similarity_NB{
	my $h_cluster_A = shift;
	my $h_cluster_B = shift;
	my $h_cluster_Mix = shift;
	my %methylTableA = calculate_methylTable_from_haps($h_cluster_A);
	my %methylTableB = calculate_methylTable_from_haps($h_cluster_B);
	my ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned)=(0,0,0);
	my ($score_A, $score_B) = (0,0);
	my ($total_hap_A, $total_hap_B) = (0,0);
	foreach my $hapString (keys(%{$h_cluster_A})){
		$total_hap_A += ${$h_cluster_A}{$hapString};
	}
	foreach my $hapString (keys(%{$h_cluster_B})){
		$total_hap_B += ${$h_cluster_B}{$hapString};
	}
	my $weight_hap_A = $total_hap_A/($total_hap_A+$total_hap_B);
	my $weight_hap_B = $total_hap_B/($total_hap_A+$total_hap_B);
	
	foreach my $hapString (keys(%{$h_cluster_Mix})){
		my $log_posterior_A=log($weight_hap_A);
		my $log_posterior_B=log($weight_hap_B);
		for(my $i=0; $i<length($hapString); $i++){
			next if(substr($hapString,$i,1) eq "N" || $methylTableA{$i}->{"depth"} eq 0  || $methylTableB{$i}->{"depth"} eq 0); 
				#don't use if not enough information from both groups
			if(substr($hapString,$i,1) =~ /[C1]/){
				my $abundance_in_A = calcSmoothAbundance($methylTableA{$i}->{"methyl"}, $methylTableA{$i}->{"depth"});
				my $abundance_in_B = calcSmoothAbundance($methylTableB{$i}->{"methyl"}, $methylTableB{$i}->{"depth"});
				$log_posterior_A += log($abundance_in_A);
				$log_posterior_B += log($abundance_in_B);
			}else{
				my $abundance_in_A = calcSmoothAbundance($methylTableA{$i}->{"depth"}-$methylTableA{$i}->{"methyl"}, $methylTableA{$i}->{"depth"});
				my $abundance_in_B = calcSmoothAbundance($methylTableB{$i}->{"depth"}-$methylTableB{$i}->{"methyl"}, $methylTableB{$i}->{"depth"});
				$log_posterior_A += log($abundance_in_A);
				$log_posterior_B += log($abundance_in_B);
			}
		}
		my $total_evidence = exp($log_posterior_A) + exp($log_posterior_B);
		if($total_evidence < $min_evidence*.1){
			$haps_unassigned += ${$h_cluster_Mix}{$hapString};
			next;
		}
		$haps_assigned_to_A += ${$h_cluster_Mix}{$hapString}*exp($log_posterior_A)/$total_evidence; # if($log_posterior_A - $log_posterior_B > 0);
		$haps_assigned_to_B += ${$h_cluster_Mix}{$hapString}*exp($log_posterior_B)/$total_evidence; # if($log_posterior_B - $log_posterior_A > 0);
		##### SCORES #####
		$score_A += $log_posterior_A - log($total_evidence);
		$score_B += $log_posterior_B - log($total_evidence);
	}
	return ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned, $score_A, $score_B);
}

#assign haplotype by similarity & abundance 
sub assign_haps_by_similarity{
	my $haps_A = shift;
	my $haps_B = shift;
	my $haps_Mix = shift;

	my ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned)=(0,0,0);
	my ($score_A, $score_B) = (0,0);

	### GENERATE HAP CLUSTER ###
	#print "$ref_hap_file_A:\n";
	my %hapClusterA = find_all_haplotype_cluster($haps_A);
	#print "$ref_hap_file_B:\n";
	my %hapClusterB = find_all_haplotype_cluster($haps_B);
	#print "$mix_hap_file:\n";
	my %hapClusterMix = find_all_haplotype_cluster($haps_Mix);
	#print "\t ", $hapClusterA{"total_hap_count"}, "\t", $hapClusterB{"total_hap_count"}, "\t", $hapClusterMix{"total_hap_count"}, "\n";
	return ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned, $score_A, $score_B) if(!$hapClusterMix{"total_hap_count"});

	my @hapClusterMix_in_A = find_similar_haplotype_binom(\%hapClusterA,\%hapClusterMix);
	my @hapClusterMix_in_B = find_similar_haplotype_binom(\%hapClusterB,\%hapClusterMix);

	my $total_hap_A = $hapClusterA{"total_hap_count"};
	my $total_hap_B = $hapClusterB{"total_hap_count"};
	my $weight_hap_A = $total_hap_A/($total_hap_A+$total_hap_B);
	my $weight_hap_B = $total_hap_B/($total_hap_A+$total_hap_B);

	foreach my $hap_info (@hapClusterMix_in_A){
		my ($hap, $score,$hapCount,$abundance_in_ref,$matched_hap) = split(/,/, $hap_info);
		$hapClusterMix{"abundance_in_ref_A"}->{$hap} = $abundance_in_ref;
		$hapClusterMix{"score_in_ref_A"}->{$hap} = $score;
		$hapClusterMix{"matched_hap_in_ref_A"}->{$hap} = $matched_hap;
		#print "RefA: $hap_info\n";
	}
	foreach my $hap_info (@hapClusterMix_in_B){
		my ($hap, $score,$hapCount,$abundance_in_ref,$matched_hap) = split(/,/, $hap_info);
		$hapClusterMix{"abundance_in_ref_B"}->{$hap} = $abundance_in_ref;
		$hapClusterMix{"score_in_ref_B"}->{$hap} = $score;
		$hapClusterMix{"matched_hap_in_ref_B"}->{$hap} = $matched_hap;
		#print "RefB: $hap_info\n";
	}
	foreach my $hap (keys(%{$hapClusterMix{"hapCounts"}})){
		$hapClusterMix{"abundance_in_ref_B"}->{$hap} ||= calcSmoothAbundance(0,$total_hap_B);
		$hapClusterMix{"abundance_in_ref_A"}->{$hap} ||= calcSmoothAbundance(0,$total_hap_A);
		$hapClusterMix{"score_in_ref_B"}->{$hap} ||= 1;
		$hapClusterMix{"score_in_ref_A"}->{$hap} ||= 1;
		$hapClusterMix{"matched_hap_in_ref_B"}->{$hap} ||= "NA";
		$hapClusterMix{"matched_hap_in_ref_A"}->{$hap} ||= "NA";
		my $abundance_in_A = $hapClusterMix{"abundance_in_ref_A"}->{$hap};
		my $abundance_in_B = $hapClusterMix{"abundance_in_ref_B"}->{$hap};

		# higher p value (scores) is worse
		if($hapClusterMix{"score_in_ref_A"}->{$hap} - $hapClusterMix{"score_in_ref_B"}->{$hap} > 0.5){
			#$abundance_in_A = calcSmoothAbundance(0,$total_hap_A); #If the difference in score is greater than 0.5, do not count for B.
		}elsif($hapClusterMix{"score_in_ref_B"}->{$hap} - $hapClusterMix{"score_in_ref_A"}->{$hap} > 0.5){
			#$abundance_in_B = calcSmoothAbundance(0,$total_hap_B); #If the difference in score is greater than 0.5, do not count for A.
		}

		my $total_evidence = $weight_hap_A*$abundance_in_A + $weight_hap_B*$abundance_in_B;
		if($total_evidence > $min_evidence){ #evidence should be at least 0.5? KZ had 0.001
			my $To_A=$hapClusterMix{"hapCounts"}->{$hap}*($weight_hap_A*$abundance_in_A/$total_evidence);
			my $To_B=$hapClusterMix{"hapCounts"}->{$hap}*($weight_hap_B*$abundance_in_B/$total_evidence);
			$haps_assigned_to_A+=$To_A;
			$haps_assigned_to_B+=$To_B;
			if($abundance_in_A != 0 && $abundance_in_B != 0){
				##### SCORES #####
				$score_A += log($weight_hap_A*$abundance_in_A/$total_evidence);
				$score_B += log($weight_hap_B*$abundance_in_B/$total_evidence);
			}
		}else{
			$haps_unassigned+=$hapClusterMix{"hapCounts"}->{$hap};
		}
	}
	return ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned, $score_A, $score_B);
}

#find haplotype groups in B that are similar to A (using binomial p value) and their relative abundance in A.
sub find_similar_haplotype_binom{
	my $h_dominantHapInfo_A = shift;
	my $h_dominantHapInfo_B = shift;
	my @similar_haps;
	my @unmatched_haps;
	my @diff_haps;	
	my @haps_A = keys(%{${$h_dominantHapInfo_A}{"haps"}});
	foreach my $hap_B (keys(%{${$h_dominantHapInfo_B}{"haps"}})){
		if(${$h_dominantHapInfo_A}{"haps"}->{$hap_B}){
			my $abundance_in_A = sprintf("%6.5f",calcSmoothAbundance(${$h_dominantHapInfo_A}{"hapCounts"}->{$hap_B},${$h_dominantHapInfo_A}{"total_hap_count"}));
			my $binomial_p = binomialCoef(0, length($hap_B))*(1-$mean_error_rate)**length($hap_B);
			push(@similar_haps, $hap_B.",".$binomial_p.",".${$h_dominantHapInfo_B}{"hapCounts"}->{$hap_B}.",$abundance_in_A,$hap_B");
		}else{			
			push(@diff_haps, $hap_B);
		}
	}
	foreach my $hap (@diff_haps){
		my ($closest_hap, $closest_distance, $max_matches)=("",100,0);		
		foreach my $hap_A (@haps_A){
			next if(length($hap) != length($hap_A));
			my $hemming_distance=0;
			my $match_bases=0;
			for(my $i=0; $i<length($hap_A); $i++){
				next if(substr($hap,$i,1) eq "N" || substr($hap_A,$i,1) eq "N");
				if(substr($hap,$i,1) ne substr($hap_A,$i,1)){
					$hemming_distance++;
				}else{
				    $match_bases++;
				}				
			}
			if($match_bases>$max_matches){
				$closest_hap = $hap_A;
				$closest_distance = $hemming_distance;
				$max_matches=$match_bases;
			}						
		}
		next if(!$closest_hap);
		my $abundance_in_A = sprintf("%6.5f", calcSmoothAbundance(${$h_dominantHapInfo_A}{"hapCounts"}->{$closest_hap},${$h_dominantHapInfo_A}{"total_hap_count"}));
		my $n = $max_matches + $closest_distance;
		my $binomial_p = 0; # the probability that there are fewer than k (closest-distance) errors
		# 1 - binomial_p is the probability that there are k or more errors.
		if(defined($binompTable{$n}->{$closest_distance-1})){
			$binomial_p = $binompTable{$n}->{$closest_distance-1};
		}else{
			my $i = 0;
			do{
				$binomial_p += binomialCoef($i, $n)*($mean_error_rate**$i)*(1-$mean_error_rate)**($n-$i);
				$i++;
			}while($i < $closest_distance);
			$binomial_p = sprintf("%.6f", $binomial_p);
			$binompTable{$n}->{$closest_distance-1} = sprintf("%.6f", $binomial_p);
		}
		if($max_matches>=4 && 1-$binomial_p >= 0.1){
			push(@similar_haps,$hap.",".$binomial_p.",".${$h_dominantHapInfo_B}{"hapCounts"}->{$hap}.",$abundance_in_A,$closest_hap");  
		}else{
			push(@unmatched_haps,"##Unmatched haps:".$hap.",".$binomial_p.",".${$h_dominantHapInfo_B}{"hapCounts"}->{$hap}.",$abundance_in_A")  
		}
	}	
	return @similar_haps;	
}

#find haplotype groups in B that are similar to A and their relative abundance in A.
sub find_similar_haplotype_z_score{
	my $h_dominantHapInfo_A = shift;
	my $h_dominantHapInfo_B = shift;
	my @similar_haps;
	my @unmatched_haps;
	my @diff_haps;	
	my @haps_A = keys(%{${$h_dominantHapInfo_A}{"haps"}});
	my $expected_distance = $mean_error_rate*length($haps_A[0]);
	foreach my $hap_B (keys(%{${$h_dominantHapInfo_B}{"haps"}})){
		if(${$h_dominantHapInfo_A}{"haps"}->{$hap_B}){
			my $abundance_in_A = sprintf("%6.5f",calcSmoothAbundance(${$h_dominantHapInfo_A}{"hapCounts"}->{$hap_B},${$h_dominantHapInfo_A}{"total_hap_count"}));
			push(@similar_haps, $hap_B.",0.001,".${$h_dominantHapInfo_B}{"hapCounts"}->{$hap_B}.",$abundance_in_A,$hap_B");
		}else{			
			push(@diff_haps, $hap_B);
		}
	}
	foreach my $hap (@diff_haps){
		my ($closest_hap, $closest_distance, $max_matches)=("",100,0);		
		foreach my $hap_A (@haps_A){
			next if(length($hap) != length($hap_A));
			my $hemming_distance=0;
			my $match_bases=0;
			for(my $i=0; $i<length($hap_A); $i++){
				next if(substr($hap,$i,1) eq "N" || substr($hap_A,$i,1) eq "N");
				if(substr($hap,$i,1) ne substr($hap_A,$i,1)){
					$hemming_distance++;
				}else{
				    $match_bases++;
				}				
			}
			if($match_bases>$max_matches){
				$closest_hap = $hap_A;
				$closest_distance = $hemming_distance;
				$max_matches=$match_bases;
			}						
		}
		next if(!$closest_hap);
		my $abundance_in_A = sprintf("%6.5f", calcSmoothAbundance(${$h_dominantHapInfo_A}{"hapCounts"}->{$closest_hap},${$h_dominantHapInfo_A}{"total_hap_count"}));
		my $Z_score = sprintf("%5.2f", abs($closest_distance-$expected_distance)/sqrt($expected_distance));
		if($Z_score<1.65 && $max_matches>=4){ #0.95
			push(@similar_haps,$hap.",".$Z_score.",".${$h_dominantHapInfo_B}{"hapCounts"}->{$hap}.",$abundance_in_A,$closest_hap");  
		}else{
			push(@unmatched_haps,"##Unmatched haps:".$hap.",".$Z_score.",".${$h_dominantHapInfo_B}{"hapCounts"}->{$hap}.",$abundance_in_A")  
		}
	}	
	#print join("\n", @unmatched_haps),"\n";
	return @similar_haps;	
}

#find all haplotype clusters
sub find_all_haplotype_cluster{
	my $h_haplotype_group=shift;
	my %dominantHapInfo;
	my $remainingHaps = 0;
	my ($dominant_hap, $haps_in_cluster, $totalHaps, $hapsCount_in_cluster);
	do{
		($dominant_hap, $haps_in_cluster, $totalHaps, $hapsCount_in_cluster)=find_dominant_haplotype_cluster($h_haplotype_group);
		if($hapsCount_in_cluster>0){
			$dominantHapInfo{"hapCounts"}->{$dominant_hap}=$hapsCount_in_cluster;
			$dominantHapInfo{"haps"}->{$dominant_hap}=$haps_in_cluster;	
			$dominantHapInfo{"total_hap_count"}+=$hapsCount_in_cluster;
			$remainingHaps = $totalHaps - $hapsCount_in_cluster; #if(!$remainingHaps);
			#$dominantHapInfo{"entropy"}-= $hapsCount_in_cluster/$totalHap*log($hapsCount_in_cluster/$totalHap)if($totalHap);
			#print "$dominant_hap\t$hapsCount_in_cluster\n";#$dominantHapInfo{"haps"}->{$dominant_hap},"\n";
			#print "##remainingHaps=$remainingHaps\n"; #hapsCount_in_cluster=$hapsCount_in_cluster\n";
		}
	}while($totalHaps > 0);
	#print "entropy=",$dominantHapInfo{"entropy"},"\n";
	return(%dominantHapInfo);
}


#find the dominant haplotype cluster
sub find_dominant_haplotype_cluster{
	my $h_haplotype_group=shift;
	my $totalHaps = 0;	
	my $dominant_hap = "";
	my $dominant_hap_count = 0;
	my $dominant_hap_valid_sites=0;
	my %methylTable;
	foreach my $hapString (keys(%{$h_haplotype_group})){
		#print "hapString:$hapString\t", ${$h_haplotype_group}{$hapString}, "\n";
		$totalHaps+=${$h_haplotype_group}{$hapString};
		my $hap_valid_sites=0;
		for(my $i=0;$i<length($hapString);$i++){
			my $c = substr($hapString,$i,1);
			$methylTable{$i}->{"depth"} = 0 if(!$methylTable{$i}->{"depth"});
			$methylTable{$i}->{"methyl"} = 0 if(!$methylTable{$i}->{"methyl"});
			$methylTable{$i}->{"depth"} += ${$h_haplotype_group}{$hapString} if($c ne "N");
			$methylTable{$i}->{"methyl"} += ${$h_haplotype_group}{$hapString} if($c =~ /[C1]/);
			$hap_valid_sites++ if($c ne "N");
			$methylTable{$i}->{"aveMethyl"} = calculate_methyl_ave_i($methylTable{$i}->{"methyl"}, $methylTable{$i}->{"depth"}) if($methylTable{$i}->{"depth"} > 0);
		}
		next if($hap_valid_sites<$dominant_hap_valid_sites);
		next if($hap_valid_sites==$dominant_hap_valid_sites && ${$h_haplotype_group}{$hapString} <=$dominant_hap_count);
		$dominant_hap = $hapString;
		$dominant_hap_count = ${$h_haplotype_group}{$hapString};
		$dominant_hap_valid_sites=$hap_valid_sites;
	}
	my $hapsCount_in_cluster=$dominant_hap_count;
	my @haps_cluster;
	push(@haps_cluster, $dominant_hap);
	my $expected_distance = $mean_error_rate*length($dominant_hap);
	foreach my $hapString (keys(%{$h_haplotype_group})){
		next if(!${$h_haplotype_group}{$hapString});
		next if($hapString eq $dominant_hap);
		my $log_p_value=0;
		my $overlapping_bases=0;
		for(my $i=0; $i<length($hapString); $i++){
			next if(!defined($methylTable{$i}->{"aveMethyl"}));
			$overlapping_bases++ if(substr($hapString,$i,1) ne "N" && substr($dominant_hap,$i,1) ne "N" );
			next if(substr($hapString,$i,1) eq "N" || 
					substr($dominant_hap,$i,1) eq "N" ||
					substr($hapString,$i,1) eq substr($dominant_hap,$i,1));
			if(substr($dominant_hap,$i,1) =~ /[C1]/){
				$log_p_value+=log(1-$methylTable{$i}->{"aveMethyl"}+0.0001);
			}else{
				$log_p_value+=log($methylTable{$i}->{"aveMethyl"}+0.0001);
			}
		}
		next if($log_p_value < -2.995732 || $overlapping_bases<3); #95% or too few overlapping bases, log(0.05) = -2.995732
		push(@haps_cluster,$hapString);
		$hapsCount_in_cluster+=${$h_haplotype_group}{$hapString};
	}
	foreach my $hapString (@haps_cluster){
		delete(${$h_haplotype_group}{$hapString});
	}
	#print "####Dominant_hap:$dominant_hap\t$hapsCount_in_cluster\t",join(",",@haps_cluster),"\n" if(!$totalHaps || !$hapsCount_in_cluster);
	return($dominant_hap, join(",",@haps_cluster), $totalHaps, $hapsCount_in_cluster);	
}

sub load_data_file{
	#1) ID
	#2) hapString
	#3) hapCount
	#4) cpgPositions
	#read ref A haplotypes
	my $file = shift;
	my %hapTable;
	my $hap_number=0;
	open(INFILE, "$file")||die("Error in opening file $file\n");
	while(my $line = <INFILE>){
		chomp($line);
		my @fields = split(/\t/, $line);
		next if(scalar(@fields)<4);
		my $probeID = $fields[0];
		my ($chr, $start, $end) = split /[:-]/, $probeID;
		my $hapString = $fields[1];
		my @bases = split "", $hapString;
		my @cpgPositions = split(/[,:]/, $fields[3]);
		next if(length($hapString)<$minHapLen || length($hapString) ne scalar(@cpgPositions));
		my $hapCount = $fields[2];
		my $start_i = 0; 
		my $end_i = length($hapString)-1;
		$hapString =~ s/C/1/g;$hapString =~ s/T/0/g; $hapString =~ s/[AG]/N/g;
		
		# remove any padded N's
		my $b = $bases[$start_i];
		while($b eq "N"){
			$b = $bases[$start_i];
			$start_i++;
		}		
		$b = $bases[$end_i];
		while($b eq "N"){
			$b = $bases[$end_i];
			$end_i--;
		}
		my $len = $end_i - $start_i + 1;
		next if($len < $minHapLen);
		$probeID = $chr.":".$cpgPositions[$start_i]."-".$cpgPositions[$end_i];
		for(my $i = $start_i; $i <= $end_i; $i++){
			$cgTable{$chr.":".$cpgPositions[$i]}->{$probeID} = 1 if($bases[$i] ne "N");
		}
		my $hapString_substr = substr($hapString, $start_i, $len);
		my @updated_cpgPositions = @cpgPositions[$start_i..$end_i];
		#fix variable length probeID
		if($hapTable{$probeID}->{"cpgPositions"}){
			my $cur_cpgPositions = $hapTable{$probeID}->{"cpgPositions"};
			# take the longest len
			if(scalar(@$cur_cpgPositions) < $len){
				# update/align all previous hapStrings
				foreach my $cur_hapString(keys %{$hapTable{$probeID}->{"haps"}}){
					my $cur_h_count = $hapTable{$probeID}->{"haps"}->{$cur_hapString};
					my $aligned_hapString = align_hapString_to_pos(\@updated_cpgPositions, $cur_cpgPositions, $cur_hapString);
					delete($hapTable{$probeID}->{"haps"}->{$cur_hapString});
					$hapTable{$probeID}->{"haps"}->{$aligned_hapString} = $cur_h_count;
				}
			}elsif(scalar(@$cur_cpgPositions) > $len){
				# update/align current hapString
				my $aligned_hapString = align_hapString_to_pos($cur_cpgPositions, \@updated_cpgPositions, $hapString_substr);
				undef(@updated_cpgPositions);
				@updated_cpgPositions = @$cur_cpgPositions;
				$hapString_substr = $aligned_hapString;
			}
		}

		$hapTable{$probeID}->{"haps"}->{$hapString_substr} += $hapCount;
		@{$hapTable{$probeID}->{"cpgPositions"}} = @updated_cpgPositions;
		$hapTable{$probeID}->{"totalHap"} += $hapCount;
		$hap_number+=$hapCount;
	}
	close(INFILE);
	return %hapTable;
	#print "$hap_number haplotypes loaded for $file\n";
}


sub binomialCoef{
        my ($k, $n) = @_;
        my $r = 1;
        do{
                $r *= $n/($n-$k);
                $n--;
        }while($n > $k);
        return $r;
}

sub calcSmoothAbundance{
	my $num_occurrences_in_docs = shift;
	my $total_docs = shift;
	return (($num_occurrences_in_docs + $smoothing)/($total_docs + 2*$smoothing));
}

sub calculate_methyl_ave_i{
	my $m = shift;
	my $d = shift;
	my $ave_m = $m/$d;
	$ave_m = $mean_error_rate if($m/$d < $mean_error_rate);
	$ave_m = 1 - $mean_error_rate if($m/$d > 1 - $mean_error_rate);
	return $ave_m;
}

#calcualte methylTable from haplotype clusters
sub calculate_methylTable_from_haps{
	my $h_cluster_group = shift;
	my %methylTable;
	foreach my $hapString (keys (%{$h_cluster_group})){
		for(my $i=0;$i<length($hapString);$i++){
			my $c = substr($hapString,$i,1);
			$methylTable{$i}->{"depth"} = 0 if(!$methylTable{$i}->{"depth"});
			$methylTable{$i}->{"methyl"} = 0 if(!$methylTable{$i}->{"methyl"});
			$methylTable{$i}->{"depth"} += ${$h_cluster_group}{$hapString} if($c ne "N");
			$methylTable{$i}->{"methyl"} += ${$h_cluster_group}{$hapString} if($c =~ /[C1]/);
		}	
	}
	return %methylTable;
}


#calculate average methylation per target
sub calc_mean_methylation{
	my $h_haplotype_group=shift;
	my %baseCounts;
	$baseCounts{"1"}=0.001;$baseCounts{"0"}=0.001;
	foreach my $hapString (keys(%{$h_haplotype_group})){
		for(my $i=0; $i<length($hapString); $i++){
			$baseCounts{substr($hapString,$i,1)}+=${$h_haplotype_group}{$hapString};
		}			
	}	
	return($baseCounts{"1"}/($baseCounts{"1"}+$baseCounts{"0"}));
}

#find all haplotype clusters
sub one_hap_per_cluster{
	my $h_haplotype_group=shift;
	my %dominantHapInfo;
	foreach my $hapString (keys(%{$h_haplotype_group})){
		$dominantHapInfo{"hapCounts"}->{$hapString}=${$h_haplotype_group}{$hapString};
		$dominantHapInfo{"haps"}->{$hapString}=$hapString;	
		$dominantHapInfo{"total_hap_count"}+=${$h_haplotype_group}{$hapString};
	}
	return%dominantHapInfo;
}


main(); 
