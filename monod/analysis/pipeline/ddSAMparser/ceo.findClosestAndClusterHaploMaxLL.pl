#!/usr/bin/perl -w
use strict;
no warnings 'recursion';

# CONSTANTS
my $smoothing = 0.001;
my $minHapLen = 4;
my $minimal_mean_mC_diff = 0.2;
my $useBinom = 1;

# optimizable parameters:
my $min_hap_counts = 50;
my $mean_error_rate = 0.03;
my $min_evidence = 0.01;
my $max_ceo_diff = -20;

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
my $use_method = $ARGV[3] ? $ARGV[3] : "abundance";

my $hap_number = 0;
my $lib_size_A = 0;
my $lib_size_B = 0;

sub main{
	
	### LOAD HAP INFO DATA FILES
	my %hapTableA = load_data_file($ref_hap_file_A);
	$lib_size_A = $hap_number;
	my %hapTableB = load_data_file($ref_hap_file_B);
	$lib_size_B = $hap_number;
	my %hapTableMix = load_data_file($mix_hap_file);
	### END LOAD HAP INFO 


	### BEGIN ASSIGNMENT
	my ($probe_used, $hap_used, $total_assigned_to_A, $total_assigned_to_B, $total_unassigned) = (0,0,0,0,0);
	my ($total_log_score_to_A, $total_log_score_to_B) = (0,0);

	open(OUT_DATA, ">MinorHaps.txt") || die("Error writing to MinorHaps.txt\n");
	foreach my $probeID (keys %hapTableMix){
		### BEGIN SEARCHING APPROPRIATE HAPLOTYPES FROM A AND B TO COMPARE WITH MIX
		my %matchedHapTable;
		my %cgOverlaps;
		$matchedHapTable{"totalHap"}->{"A"} = 0;
		$matchedHapTable{"totalHap"}->{"B"} = 0;
		my ($chr,$chr_start,$chr_end) = split /[:-]/, $probeID;
		my @cpgPositions = @{$hapTableMix{$probeID}->{"cpgPositions"}};
		next if(scalar(@cpgPositions) < $minHapLen || $chr =~ /[XYM]/);
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
		### END SEARCH

		### DECIDE IF PROBEID/REGION IS APPROPRIATE
		next if( $matchedHapTable{"totalHap"}->{"A"} < $min_hap_counts || $matchedHapTable{"totalHap"}->{"B"} < $min_hap_counts);
		my $refA_mean_methylation = sprintf("%4.3f", calc_mean_methylation(\%{$matchedHapTable{"hapCounts"}->{"A"}}));
		my $refB_mean_methylation = sprintf("%4.3f", calc_mean_methylation(\%{$matchedHapTable{"hapCounts"}->{"B"}}));
		next if(abs($refB_mean_methylation-$refA_mean_methylation)<$minimal_mean_mC_diff);
		next if(calc_ceo(\%{$matchedHapTable{"hapCounts"}->{"A"}}, \%{$matchedHapTable{"hapCounts"}->{"B"}}) > $max_ceo_diff);
		$probe_used++;
		$hap_used += $hapTableMix{$probeID}->{"totalHap"};
		$read_counts_by_chr{$chr}->{"refA"} = 0 if(!$read_counts_by_chr{$chr}->{"refA"});
		$read_counts_by_chr{$chr}->{"refB"} = 0 if(!$read_counts_by_chr{$chr}->{"refB"});

		### ASSIGN HAPLOTYPES 
		my ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned, $score_A, $score_B) = assign_switch(\%{$matchedHapTable{"hapCounts"}->{"A"}}, \%{$matchedHapTable{"hapCounts"}->{"B"}}, \%{$hapTableMix{$probeID}->{"haps"}});
		print OUT_DATA "$probeID\t$haps_assigned_to_B\t$haps_unassigned\t", $hapTableMix{$probeID}->{"totalHap"}, "\n";
		$read_counts_by_chr{$chr}->{"refA"}+=$haps_assigned_to_A;
		$read_counts_by_chr{$chr}->{"refB"}+=$haps_assigned_to_B;
		$read_counts_by_chr{$chr}->{"unassigned"}+=$haps_unassigned;
		$total_log_score_to_A += $score_A;
		$total_log_score_to_B += $score_B;
		undef(%matchedHapTable);
		undef(%cgOverlaps);
	}
	close(OUT_DATA);
	### END ASSIGNMENT

	# BEGIN OUTPUT
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
	# END OUTPUT

	undef(%hapTableA);
	undef(%hapTableB);
	undef(%hapTableMix);
	undef(%read_counts_by_chr);
}

sub assign_switch{
	my $h_haps_A = shift;
	my $h_haps_B = shift;
	my $h_haps_Mix = shift;
	if($use_method eq "abundance"){
		return assign_haps_by_abundance($h_haps_A, $h_haps_B, $h_haps_Mix);
	}elsif($use_method eq "bernoulliNB_cg"){
		return assign_haps_by_bernoulliNB_cg($h_haps_A, $h_haps_B, $h_haps_Mix);
	}elsif($use_method eq "bernoulliNB"){
		return assign_haps_by_bernoulliNB($h_haps_A, $h_haps_B, $h_haps_Mix);
	}elsif($use_method eq "multinomialNB"){
		return assign_haps_by_multinomialNB($h_haps_A, $h_haps_B, $h_haps_Mix);
	}else{
		return assign_haps_by_abundance($h_haps_A, $h_haps_B, $h_haps_Mix);
	}
}

#assign haplotype by NB using CG sites as word
sub assign_haps_by_bernoulliNB_cg{
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

	foreach my $hapString (keys(%{$h_cluster_Mix})){
		my $log_posterior_A=log($total_hap_A/($total_hap_A+$total_hap_B));
		my $log_posterior_B=log($total_hap_B/($total_hap_A+$total_hap_B));
		for(my $i=0; $i<length($hapString); $i++){
			next if(substr($hapString,$i,1) eq "N" || $methylTableA{$i}->{"depth"} eq 0  || $methylTableB{$i}->{"depth"} eq 0); 
				#don't use if not enough information from both groups
			if(substr($hapString,$i,1) =~ /[C1]/){
				my $abundance_in_A = calcSmoothAbundance($methylTableA{$i}->{"methyl"}, $methylTableA{$i}->{"depth"},2);
				my $abundance_in_B = calcSmoothAbundance($methylTableB{$i}->{"methyl"}, $methylTableB{$i}->{"depth"},2);
				$log_posterior_A += log($abundance_in_A);
				$log_posterior_B += log($abundance_in_B);
			}else{
				my $abundance_in_A = calcSmoothAbundance($methylTableA{$i}->{"depth"}-$methylTableA{$i}->{"methyl"}, $methylTableA{$i}->{"depth"},2);
				my $abundance_in_B = calcSmoothAbundance($methylTableB{$i}->{"depth"}-$methylTableB{$i}->{"methyl"}, $methylTableB{$i}->{"depth"},2);
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

#assign haplotype by similarity assign Naive Bayes scores using multinomial NB
sub assign_haps_by_multinomialNB{
	my $h_haps_A = shift;
	my $h_haps_B = shift;
	my $h_haps_Mix = shift;
	my ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned)=(0,0,0);
	my ($score_A, $score_B) = (0,0);

	my @haps_A = keys %{$h_haps_A};
	my @haps_B = keys %{$h_haps_B};
	my @haps_Mix = keys %{$h_haps_Mix};
	my %haps_Bag;

	my $mix_hap_count = 0;
	### GENERATE BAG OF HAPS CLUSTER
	foreach my $hapString (@haps_A){
		$haps_Bag{$hapString} += ${$h_haps_A}{$hapString};
	}
	foreach my $hapString (@haps_B){
		$haps_Bag{$hapString} += ${$h_haps_B}{$hapString};
	}
	foreach my $hapString (@haps_Mix){
		$haps_Bag{$hapString} += ${$h_haps_Mix}{$hapString};
		$mix_hap_count += ${$h_haps_Mix}{$hapString};
	}
	my %hapClusterBag = find_all_haplotype_cluster(\%haps_Bag);
	undef %haps_Bag;
	return ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned, $score_A, $score_B) if(!$hapClusterBag{"total_hap_count"});

	### FOR EACH HAP in A, B, and Mix, find the best match in Haps Bag
	my @hap_features = keys %{$hapClusterBag{"hapCounts"}};
	my $total_hap_words = $hapClusterBag{"total_hap_count"};
	my $total_diff_haps = scalar(@hap_features);
	$hapClusterBag{"total_hap_in_A"} = 0;
	$hapClusterBag{"total_hap_in_B"} = 0;
	$hapClusterBag{"total_hap_in_Mix"} = 0;
	my @hapClusterA_in_Bag;
	my @hapClusterB_in_Bag;
	my @hapClusterMix_in_Bag;
	my @diff_haps_A;
	my @diff_haps_B;
	my @diff_haps_Mix;
	foreach my $hap_word (@haps_A){
		push(@hapClusterA_in_Bag, "0.001,".$hap_word.",".$hap_word) if($hapClusterBag{"hapCounts"}->{$hap_word});
		push(@diff_haps_A, $hap_word) if(!$hapClusterBag{"hapCounts"}->{$hap_word});
	}
	foreach my $hap_word (@haps_B){
		push(@hapClusterB_in_Bag, "0.001,".$hap_word.",".$hap_word) if($hapClusterBag{"hapCounts"}->{$hap_word});
		push(@diff_haps_B, $hap_word) if(!$hapClusterBag{"hapCounts"}->{$hap_word});
	}
	foreach my $hap_word (@haps_Mix){
		push(@hapClusterMix_in_Bag, "0.001,".$hap_word.",".$hap_word) if($hapClusterBag{"hapCounts"}->{$hap_word});
		push(@diff_haps_Mix, $hap_word) if(!$hapClusterBag{"hapCounts"}->{$hap_word});
	}
	foreach my $hap_info (@hapClusterA_in_Bag, find_similar_haplotype(\@hap_features, \@diff_haps_A)){
		my ($score,$hap,$matched_hap) = split(/,/, $hap_info);
		print "\t A : $hap_info\n" if(!$hapClusterBag{"hapCounts"}->{$matched_hap});
		$hapClusterBag{"count_in_A"}->{$matched_hap} += ${$h_haps_A}{$hap};
		$hapClusterBag{"total_hap_in_A"} += ${$h_haps_A}{$hap};
	}
	foreach my $hap_info (@hapClusterB_in_Bag, find_similar_haplotype(\@hap_features, \@diff_haps_B)){
		my ($score,$hap,$matched_hap) = split(/,/, $hap_info);
		print "\t B : $hap_info\n" if(!$hapClusterBag{"hapCounts"}->{$matched_hap});
		$hapClusterBag{"count_in_B"}->{$matched_hap} += ${$h_haps_B}{$hap};
		$hapClusterBag{"total_hap_in_B"} += ${$h_haps_B}{$hap};
	}
	foreach my $hap_info (@hapClusterMix_in_Bag, find_similar_haplotype(\@hap_features, \@diff_haps_Mix)){
		my ($score,$hap,$matched_hap) = split(/,/, $hap_info);
		print "\t Mix : $hap_info\n" if(!$hapClusterBag{"hapCounts"}->{$matched_hap});
		$hapClusterBag{"count_in_Mix"}->{$matched_hap} += ${$h_haps_Mix}{$hap};
		$hapClusterBag{"total_hap_in_Mix"} += ${$h_haps_Mix}{$hap};
	}
	my $total_hap_A = $hapClusterBag{"total_hap_in_A"};
	my $total_hap_B = $hapClusterBag{"total_hap_in_B"};
	$haps_unassigned += $mix_hap_count - $hapClusterBag{"total_hap_in_Mix"};
	return ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned, $score_A, $score_B) 
		if($hapClusterBag{"total_hap_in_Mix"} == 0 || $total_hap_A < $min_hap_counts || $total_hap_B < $min_hap_counts);
	my $log_posterior_A = log($total_hap_A/($total_hap_A+$total_hap_B));
	my $log_posterior_B = log($total_hap_B/($total_hap_A+$total_hap_B));
	my $num_features_in_mix = 0;
	foreach my $hap_word (@hap_features){
		$hapClusterBag{"count_in_A"}->{$hap_word} ||= 0;
		$hapClusterBag{"count_in_B"}->{$hap_word} ||= 0;
		$hapClusterBag{"count_in_Mix"}->{$hap_word} ||= 0;
		my $theta_A = calcSmoothAbundance($hapClusterBag{"count_in_A"}->{$hap_word}, $total_hap_A, $total_diff_haps);
		my $theta_B = calcSmoothAbundance($hapClusterBag{"count_in_B"}->{$hap_word}, $total_hap_B, $total_diff_haps);
		if($theta_A + $theta_B < $min_evidence){
			$haps_unassigned += $hapClusterBag{"count_in_Mix"}->{$hap_word};
		}else{
			$haps_assigned_to_A += $hapClusterBag{"count_in_Mix"}->{$hap_word}*$theta_A/($theta_A+$theta_B);
			$haps_assigned_to_B += $hapClusterBag{"count_in_Mix"}->{$hap_word}*$theta_B/($theta_A+$theta_B);
		}
		if($hapClusterBag{"count_in_Mix"}->{$hap_word} > 0 && $hapClusterBag{"count_in_Mix"}->{$hap_word} < 1000000){
			$log_posterior_A += $hapClusterBag{"count_in_Mix"}->{$hap_word} * log($theta_A) - logfactorial($hapClusterBag{"count_in_Mix"}->{$hap_word});
			$log_posterior_B += $hapClusterBag{"count_in_Mix"}->{$hap_word} * log($theta_B) - logfactorial($hapClusterBag{"count_in_Mix"}->{$hap_word});
		}
	}
	return ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned, $score_A, $score_B) if($hapClusterBag{"total_hap_in_Mix"} > 1000000);
	$log_posterior_A += logfactorial($hapClusterBag{"total_hap_in_Mix"});
	$log_posterior_B += logfactorial($hapClusterBag{"total_hap_in_Mix"});

	my $total_evidence = exp($log_posterior_A) + exp($log_posterior_B);
	##### SCORES #####
	$score_A += $log_posterior_A - log($total_evidence);
	$score_B += $log_posterior_B - log($total_evidence);
	undef %hapClusterBag;
	return ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned, $score_A, $score_B);
}

#assign haplotype by similarity assign Naive Bayes scores using bernoulliNB
sub assign_haps_by_bernoulliNB{
	my $h_haps_A = shift;
	my $h_haps_B = shift;
	my $h_haps_Mix = shift;
	my ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned)=(0,0,0);
	my ($score_A, $score_B) = (0,0);

	my @haps_A = keys %{$h_haps_A};
	my @haps_B = keys %{$h_haps_B};
	my @haps_Mix = keys %{$h_haps_Mix};
	my %haps_Bag;

	my $mix_hap_count = 0;
	### GENERATE BAG OF HAPS CLUSTER
	foreach my $hapString (@haps_A){
		$haps_Bag{$hapString} += ${$h_haps_A}{$hapString};
	}
	foreach my $hapString (@haps_B){
		$haps_Bag{$hapString} += ${$h_haps_B}{$hapString};
	}
	foreach my $hapString (@haps_Mix){
		$haps_Bag{$hapString} += ${$h_haps_Mix}{$hapString};
		$mix_hap_count += ${$h_haps_Mix}{$hapString};
	}
	my %hapClusterBag = find_all_haplotype_cluster(\%haps_Bag);
	undef %haps_Bag;
	return ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned, $score_A, $score_B) if(!$hapClusterBag{"total_hap_count"});

	### FOR EACH HAP in A, B, and Mix, find the best match in Haps Bag
	my @hap_features = keys %{$hapClusterBag{"hapCounts"}};
	my $total_hap_words = $hapClusterBag{"total_hap_count"};
	my $total_diff_haps = scalar(@hap_features);
	$hapClusterBag{"total_hap_in_A"} = 0;
	$hapClusterBag{"total_hap_in_B"} = 0;
	$hapClusterBag{"total_hap_in_Mix"} = 0;
	my @hapClusterA_in_Bag;
	my @hapClusterB_in_Bag;
	my @hapClusterMix_in_Bag;
	my @diff_haps_A;
	my @diff_haps_B;
	my @diff_haps_Mix;
	foreach my $hap_word (@haps_A){
		push(@hapClusterA_in_Bag, "0.001,".$hap_word.",".$hap_word) if($hapClusterBag{"hapCounts"}->{$hap_word});
		push(@diff_haps_A, $hap_word) if(!$hapClusterBag{"hapCounts"}->{$hap_word});
	}
	foreach my $hap_word (@haps_B){
		push(@hapClusterB_in_Bag, "0.001,".$hap_word.",".$hap_word) if($hapClusterBag{"hapCounts"}->{$hap_word});
		push(@diff_haps_B, $hap_word) if(!$hapClusterBag{"hapCounts"}->{$hap_word});
	}
	foreach my $hap_word (@haps_Mix){
		push(@hapClusterMix_in_Bag, "0.001,".$hap_word.",".$hap_word) if($hapClusterBag{"hapCounts"}->{$hap_word});
		push(@diff_haps_Mix, $hap_word) if(!$hapClusterBag{"hapCounts"}->{$hap_word});
	}
	foreach my $hap_info (@hapClusterA_in_Bag, find_similar_haplotype(\@hap_features, \@diff_haps_A)){
		my ($score,$hap,$matched_hap) = split(/,/, $hap_info);
		print "\t A : $hap_info\n" if(!$hapClusterBag{"hapCounts"}->{$matched_hap});
		$hapClusterBag{"count_in_A"}->{$matched_hap} += ${$h_haps_A}{$hap};
		$hapClusterBag{"total_hap_in_A"} += ${$h_haps_A}{$hap};
	}
	foreach my $hap_info (@hapClusterB_in_Bag, find_similar_haplotype(\@hap_features, \@diff_haps_B)){
		my ($score,$hap,$matched_hap) = split(/,/, $hap_info);
		print "\t B : $hap_info\n" if(!$hapClusterBag{"hapCounts"}->{$matched_hap});
		$hapClusterBag{"count_in_B"}->{$matched_hap} += ${$h_haps_B}{$hap};
		$hapClusterBag{"total_hap_in_B"} += ${$h_haps_B}{$hap};
	}
	foreach my $hap_info (@hapClusterMix_in_Bag, find_similar_haplotype(\@hap_features, \@diff_haps_Mix)){
		my ($score,$hap,$matched_hap) = split(/,/, $hap_info);
		print "\t Mix : $hap_info\n" if(!$hapClusterBag{"hapCounts"}->{$matched_hap});
		$hapClusterBag{"count_in_Mix"}->{$matched_hap} += ${$h_haps_Mix}{$hap};
		$hapClusterBag{"total_hap_in_Mix"} += ${$h_haps_Mix}{$hap};
	}
	my $total_hap_A = $hapClusterBag{"total_hap_in_A"};
	my $total_hap_B = $hapClusterBag{"total_hap_in_B"};
	$haps_unassigned += $mix_hap_count - $hapClusterBag{"total_hap_in_Mix"};
	return ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned, $score_A, $score_B) 
		if($hapClusterBag{"total_hap_in_Mix"} == 0 || $total_hap_A < $min_hap_counts || $total_hap_B < $min_hap_counts);
	my $log_posterior_A = log($total_hap_A/($total_hap_A+$total_hap_B));
	my $log_posterior_B = log($total_hap_B/($total_hap_A+$total_hap_B));
	foreach my $hap_word (@hap_features){
		$hapClusterBag{"count_in_A"}->{$hap_word} ||= 0;
		$hapClusterBag{"count_in_B"}->{$hap_word} ||= 0;
		$hapClusterBag{"count_in_Mix"}->{$hap_word} ||= 0;
		my $theta_A = calcSmoothAbundance($hapClusterBag{"count_in_A"}->{$hap_word}, $total_hap_A, $total_diff_haps);
		my $theta_B = calcSmoothAbundance($hapClusterBag{"count_in_B"}->{$hap_word}, $total_hap_B, $total_diff_haps);
		if($theta_A + $theta_B < $min_evidence){
			$haps_unassigned += $hapClusterBag{"count_in_Mix"}->{$hap_word};
		}else{
			$haps_assigned_to_A += $hapClusterBag{"count_in_Mix"}->{$hap_word} * $theta_A/($theta_A+$theta_B);
			$haps_assigned_to_B += $hapClusterBag{"count_in_Mix"}->{$hap_word} * $theta_B/($theta_A+$theta_B);
		}
		if($hapClusterBag{"count_in_Mix"}->{$hap_word} > 0){
			$log_posterior_A += log($theta_A);
			$log_posterior_B += log($theta_B);
		}else{
			$log_posterior_A += log(1 - $theta_A);
			$log_posterior_B += log(1 - $theta_B);
		}
	}
	my $total_evidence = exp($log_posterior_A) + exp($log_posterior_B);
	##### SCORES #####
	$score_A += $log_posterior_A - log($total_evidence);
	$score_B += $log_posterior_B - log($total_evidence);
	undef %hapClusterBag;
	return ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned, $score_A, $score_B);
}

#assign haplotype by abundance 
sub assign_haps_by_abundance{
	my $h_haps_A = shift;
	my $h_haps_B = shift;
	my $h_haps_Mix = shift;
	
	my ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned)=(0,0,0);
	my ($score_A, $score_B) = (0,0);

	### GENERATE HAP CLUSTER 
	#print "$ref_hap_file_A:\n";
	my %hapClusterA = find_all_haplotype_cluster($h_haps_A);
	#print "$ref_hap_file_B:\n";
	my %hapClusterB = find_all_haplotype_cluster($h_haps_B);
	#print "$mix_hap_file:\n";
	my %hapClusterMix = find_all_haplotype_cluster($h_haps_Mix);
	#print "\t ", $hapClusterA{"total_hap_count"}, "\t", $hapClusterB{"total_hap_count"}, "\t", $hapClusterMix{"total_hap_count"}, "\n";
	return ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned, $score_A, $score_B) if(!$hapClusterMix{"total_hap_count"});

	my $total_hap_A = $hapClusterA{"total_hap_count"};
	my $total_hap_B = $hapClusterB{"total_hap_count"};
	my $weight_hap_A = $total_hap_A/($total_hap_A+$total_hap_B);
	my $weight_hap_B = $total_hap_B/($total_hap_A+$total_hap_B);

	my @hapClusterMix_in_A;
	my @hapClusterMix_in_B;
	my @haps_A = keys(%{$hapClusterA{"hapCounts"}});
	my @haps_B = keys(%{$hapClusterB{"hapCounts"}});
	my @diff_haps_A;
	my @diff_haps_B;
	foreach my $hap_Mix (keys %{$hapClusterMix{"hapCounts"}}){
		push(@hapClusterMix_in_A, "0.001,".$hap_Mix.",".$hap_Mix) if($hapClusterA{"hapCounts"}->{$hap_Mix});
		push(@hapClusterMix_in_B, "0.001,".$hap_Mix.",".$hap_Mix) if($hapClusterB{"hapCounts"}->{$hap_Mix});
		push(@diff_haps_A, $hap_Mix) if(!$hapClusterA{"hapCounts"}->{$hap_Mix});
		push(@diff_haps_B, $hap_Mix) if(!$hapClusterB{"hapCounts"}->{$hap_Mix});
	}
	foreach my $hap_info (@hapClusterMix_in_A, find_similar_haplotype(\@haps_A, \@diff_haps_A)){
		my ($score,$hap,$matched_hap) = split(/,/, $hap_info);
		#print "\t A : $hap_info\n" if(!$hapClusterA{"hapCounts"}->{$matched_hap});
		$hapClusterMix{"abundance_in_ref_A"}->{$hap} = calcSmoothAbundance($hapClusterA{"hapCounts"}->{$matched_hap},$total_hap_A,2);
		$hapClusterMix{"score_in_ref_A"}->{$hap} = $score;
		$hapClusterMix{"matched_hap_in_ref_A"}->{$hap} = $matched_hap;
	}
	foreach my $hap_info (@hapClusterMix_in_B, find_similar_haplotype(\@haps_B, \@diff_haps_B)){
		my ($score,$hap,$matched_hap) = split(/,/, $hap_info);
		#print "\t B : $hap_info\n" if(!$hapClusterB{"hapCounts"}->{$matched_hap});
		$hapClusterMix{"abundance_in_ref_B"}->{$hap} = calcSmoothAbundance($hapClusterB{"hapCounts"}->{$matched_hap},$total_hap_B,2);
		$hapClusterMix{"score_in_ref_B"}->{$hap} = $score;
		$hapClusterMix{"matched_hap_in_ref_B"}->{$hap} = $matched_hap;
	}

	foreach my $hap (keys(%{$hapClusterMix{"hapCounts"}})){
		$hapClusterMix{"abundance_in_ref_B"}->{$hap} ||= calcSmoothAbundance(0,$total_hap_B,2);
		$hapClusterMix{"abundance_in_ref_A"}->{$hap} ||= calcSmoothAbundance(0,$total_hap_A,2);
		$hapClusterMix{"matched_hap_in_ref_B"}->{$hap} ||= "NA";
		$hapClusterMix{"matched_hap_in_ref_A"}->{$hap} ||= "NA";
		my $abundance_in_A = $hapClusterMix{"abundance_in_ref_A"}->{$hap};
		my $abundance_in_B = $hapClusterMix{"abundance_in_ref_B"}->{$hap};
		if($useBinom){
			$hapClusterMix{"score_in_ref_B"}->{$hap} ||= 1;
			$hapClusterMix{"score_in_ref_A"}->{$hap} ||= 1;
		}else{
			$hapClusterMix{"score_in_ref_B"}->{$hap} ||= 10;
			$hapClusterMix{"score_in_ref_A"}->{$hap} ||= 10;
		}
		# higher p value (scores) is worse
		#if($hapClusterMix{"score_in_ref_A"}->{$hap} - $hapClusterMix{"score_in_ref_B"}->{$hap} > 0.5){
			#$abundance_in_A = calcSmoothAbundance(0,$total_hap_A,2); #If the difference in score is greater than 0.5, do not count for B.
		#}elsif($hapClusterMix{"score_in_ref_B"}->{$hap} - $hapClusterMix{"score_in_ref_A"}->{$hap} > 0.5){
			#$abundance_in_B = calcSmoothAbundance(0,$total_hap_B,2); #If the difference in score is greater than 0.5, do not count for A.
		#}
		my $total_evidence = $weight_hap_A*$abundance_in_A + $weight_hap_B*$abundance_in_B;
		if($total_evidence > $min_evidence){ #evidence should be at least 0.5? KZ had 0.001
			my $To_A=$hapClusterMix{"hapCounts"}->{$hap}*($weight_hap_A * $abundance_in_A/$total_evidence);
			my $To_B=$hapClusterMix{"hapCounts"}->{$hap}*($weight_hap_B * $abundance_in_B/$total_evidence);
			$haps_assigned_to_A+=$To_A;
			$haps_assigned_to_B+=$To_B;
			if($abundance_in_A != 0 && $abundance_in_B != 0){
				##### SCORES #####
				$score_A += log($weight_hap_A * $abundance_in_A/$total_evidence);
				$score_B += log($weight_hap_B * $abundance_in_B/$total_evidence);
			}
		}else{
			$haps_unassigned+=$hapClusterMix{"hapCounts"}->{$hap};
		}
	}
	return ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned, $score_A, $score_B);
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

#find haplotype groups in B that are similar to A and their relative abundance in A.
sub find_similar_haplotype{
	my $haps_A = shift;
	my $haps_B = shift;
	my @similar_haps;
	my @unmatched_haps;
	foreach my $hap (@{$haps_B}){
		my ($closest_hap, $closest_distance, $max_matches)=("",100,0);		
		foreach my $hap_A (@{$haps_A}){
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
		my $n = $max_matches + $closest_distance;
		if(!$useBinom){
			my $expected_distance = $mean_error_rate*length($hap);
			my $Z_score = sprintf("%5.2f", abs($closest_distance-$expected_distance)/sqrt($expected_distance));
			if($Z_score<1.65 && $max_matches>=4){ #0.95
				push(@similar_haps, $Z_score.",".$hap.",".$closest_hap); 
			}else{
				push(@unmatched_haps,"##Unmatched haps:".$Z_score.",".$hap.",".$closest_hap);
			}
		}else{
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
				push(@similar_haps, $binomial_p.",".$hap.",".$closest_hap);  
			}else{
				push(@unmatched_haps,"##Unmatched haps:".$binomial_p.",".$hap.",".$closest_hap);
			}
		}
	}	
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
			next if($methylTable{$i}->{"depth"} == 0);
			my $ave_methyl = calculate_methyl_ave_i($methylTable{$i}->{"methyl"}, $methylTable{$i}->{"depth"});
			$overlapping_bases++ if(substr($hapString,$i,1) ne "N" && substr($dominant_hap,$i,1) ne "N" );
			next if(substr($hapString,$i,1) eq "N" || 
					substr($dominant_hap,$i,1) eq "N" ||
					substr($hapString,$i,1) eq substr($dominant_hap,$i,1));
			if(substr($dominant_hap,$i,1) =~ /[C1]/){
				$log_p_value+=log(1-$ave_methyl+0.0001);
			}else{
				$log_p_value+=log($ave_methyl+0.0001);
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
	$hap_number=0;
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
	my $total_features = shift;
	return (($num_occurrences_in_docs + $smoothing)/($total_docs + $total_features*$smoothing));
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

#calculate the combinatorial entropy optimization for locus
sub calc_ceo{
	my $h_haplotype_group_A=shift;
	my $h_haplotype_group_B=shift;
	my %allelesTable;
	my $N_A = 0;
	my $N_B = 0;
	foreach my $hapString (keys(%{$h_haplotype_group_A})){
		$allelesTable{$hapString}->{"A"} = ${$h_haplotype_group_A}{$hapString};
		$N_A+=${$h_haplotype_group_A}{$hapString};
	}
	foreach my $hapString (keys(%{$h_haplotype_group_B})){
		$allelesTable{$hapString}->{"B"} = ${$h_haplotype_group_B}{$hapString};
		$N_B+=${$h_haplotype_group_B}{$hapString};
	}
	my $normalized_N_A = 200*($N_A/$lib_size_A)/($N_A/$lib_size_A+$N_B/$lib_size_B);
	my $normalized_N_B = 200*($N_B/$lib_size_B)/($N_A/$lib_size_B+$N_B/$lib_size_B);
	my $foreground_subtract = 0;
	my $background_subtract = 0;
	foreach my $hapString (keys %allelesTable){
		my ($a, $b) = (0,0);
		$a = $allelesTable{$hapString}->{"A"} if($allelesTable{$hapString}->{"A"});
		$b = $allelesTable{$hapString}->{"B"} if($allelesTable{$hapString}->{"B"});
		$a = $a*$normalized_N_A/$N_A;
		$b = $b*$normalized_N_B/$N_B;
		my $N_i_j = $a + $b;
		$foreground_subtract -= logfactorial_gamma($a);
		$foreground_subtract -= logfactorial_gamma($b);
		$background_subtract -= logfactorial_gamma($normalized_N_A*$N_i_j/200);
		$background_subtract -= logfactorial_gamma($normalized_N_B*$N_i_j/200);
	}
	return $foreground_subtract - $background_subtract;
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

sub logfactorial{
	my $n = shift;
	return 0 if($n <= 1);
	return log($n)+logfactorial($n-1);
}

sub logfactorial_gamma{
	my $n = shift;
	my $z = $n + 1;
	# only works for positive real numbers.
	my @q = (75122.6331530, 80916.6278952, 36308.2951477, 8687.24529705, 1168.92649479, 83.8676043424, 2.50662827511);
	my $log_factorial = ($z+0.5)*log($z+5.5) - $z - 5.5;
	my $denom_sum = 0;
	for(my $i = 0; $i < 7; $i++){
		$log_factorial -= log($z+$i);
		$denom_sum += $q[$i]*$z^$i;
	}
	$log_factorial += log($denom_sum);
	return $log_factorial;
}

main(); 
