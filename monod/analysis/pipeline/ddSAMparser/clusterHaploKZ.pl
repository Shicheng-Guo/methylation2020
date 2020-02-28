#!/usr/bin/perl -w
use strict;

use Statistics::Regression;

my $mean_error_rate = 0.1;
my $minHapLen = 4;
my %hapTable;

# probeID=>hapCounts=>sampleID
#	 =>CpgPositions
#		 =>hapCounts
#		 =>totalHap

my $ref_hap_file_A = $ARGV[0];
my $ref_hap_file_B = $ARGV[1];
my $mix_hap_file = $ARGV[2];
my $minimal_mean_mC_diff = 0.2;
my $min_hap_counts = $ARGV[3];
my %read_counts_by_chr;

sub main{
	load_data_files();
	my ($probe_used, $hap_used, $total_assigned_to_A, $total_assigned_to_B, $total_unassigned) = (0,0,0,0,0);
	foreach my $probeID (keys(%hapTable)){
		next if($probeID =~ /chr[XYM]/);
		next if(!$hapTable{$probeID}->{"cpgPositions"});
		my @cpgPositions = @{$hapTable{$probeID}->{"cpgPositions"}};		
		next if(scalar(@cpgPositions)<$minHapLen || 
			!$hapTable{$probeID}->{"hapCounts"}->{"A"} ||
			!$hapTable{$probeID}->{"hapCounts"}->{"B"} ||
			!$hapTable{$probeID}->{"hapCounts"}->{"mix"} ||
			$hapTable{$probeID}->{"totalHap"}->{"A"}<$min_hap_counts || 
			$hapTable{$probeID}->{"totalHap"}->{"B"}<$min_hap_counts || 
			$hapTable{$probeID}->{"totalHap"}->{"mix"}<$min_hap_counts);
		my $refA_mean_methylation = sprintf("%4.3f", calc_mean_methylation(\%{$hapTable{$probeID}->{"hapCounts"}->{"A"}}));
		my $refB_mean_methylation = sprintf("%4.3f", calc_mean_methylation(\%{$hapTable{$probeID}->{"hapCounts"}->{"B"}}));
		next if(abs($refB_mean_methylation-$refA_mean_methylation)<$minimal_mean_mC_diff);
		$probe_used++;
		#print "$probeID\t$ref_hap_file_A mean mC: $refA_mean_methylation\t$ref_hap_file_B mean mC: $refB_mean_methylation\n";
		my ($chr, $chr_start, $chr_end) = split (/[\-:]/, $probeID);
		$hap_used+= $hapTable{$probeID}->{"totalHap"}->{"mix"};
		$read_counts_by_chr{$chr}->{"refA"} = 0 if(!$read_counts_by_chr{$chr}->{"refA"});
		$read_counts_by_chr{$chr}->{"refB"} = 0 if(!$read_counts_by_chr{$chr}->{"refB"});
		#print "$ref_hap_file_A:\n";
		my %hapClusterA = find_all_haplotype_cluster(\%{$hapTable{$probeID}->{"hapCounts"}->{"A"}});
		#print "$ref_hap_file_B:\n";
		my %hapClusterB = find_all_haplotype_cluster(\%{$hapTable{$probeID}->{"hapCounts"}->{"B"}});
		#print "$mix_hap_file:\n";
		my %hapClusterMix = find_all_haplotype_cluster(\%{$hapTable{$probeID}->{"hapCounts"}->{"mix"}});
		#print "\t ", $hapClusterA{"total_hap_count"}, "\t", $hapClusterB{"total_hap_count"}, "\t", $hapClusterMix{"total_hap_count"}, "\n";
		next if(!$hapClusterMix{"total_hap_count"});
		my ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned) = assign_haps_by_similarity(\%hapClusterA, \%hapClusterB, \%hapClusterMix);
		$read_counts_by_chr{$chr}->{"refA"}+=$haps_assigned_to_A;
		$read_counts_by_chr{$chr}->{"refB"}+=$haps_assigned_to_B;
		$read_counts_by_chr{$chr}->{"unassigned"}+=$haps_unassigned;
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
	$total_assigned_to_A = sprintf("%.0f", $total_assigned_to_A);
	$total_assigned_to_B = sprintf("%.0f", $total_assigned_to_B);
	print "RefA=$total_assigned_to_A RefB=$total_assigned_to_B Unassigned=$total_unassigned ";
	print "Percent $ref_hap_file_A=$percent_A (%) Percent $ref_hap_file_B=$percent_B (%)\n";
}

#assign haplotype by similarity & abundance 
sub assign_haps_by_similarity(){
	my $h_cluster_A = shift;
	my $h_cluster_B = shift;
	my $h_cluster_Mix = shift;
	my $reg_hap = Statistics::Regression->new(" Hap freq regression ", [ "const", "refA", "refB" ] );
	my @hapClusterMix_in_A = find_similar_haplotype($h_cluster_A,$h_cluster_Mix);
	my @hapClusterMix_in_B = find_similar_haplotype($h_cluster_B,$h_cluster_Mix);
	foreach my $hap_info (@hapClusterMix_in_A){
		my ($hap, $z_score,$hapCount,$abundance_in_ref,$matched_hap) = split(/,/, $hap_info);
		${$h_cluster_Mix}{"abundance_in_ref_A"}->{$hap} = $abundance_in_ref;
		${$h_cluster_Mix}{"z_score_in_ref_A"}->{$hap} = $z_score;
		${$h_cluster_Mix}{"matched_hap_in_ref_A"}->{$hap} = $matched_hap;
		#print "RefA: $hap_info\n";
	}
	foreach my $hap_info (@hapClusterMix_in_B){
		my ($hap, $z_score,$hapCount,$abundance_in_ref,$matched_hap) = split(/,/, $hap_info);
		${$h_cluster_Mix}{"abundance_in_ref_B"}->{$hap} = $abundance_in_ref;
		${$h_cluster_Mix}{"z_score_in_ref_B"}->{$hap} = $z_score;
		${$h_cluster_Mix}{"matched_hap_in_ref_B"}->{$hap} = $matched_hap;
		#print "RefB: $hap_info\n";
	}
	my ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned)=(0,0,0);
	foreach my $hap (keys(%{${$h_cluster_Mix}{"hapCounts"}})){
		${$h_cluster_Mix}{"abundance_in_ref_B"}->{$hap} ||= 0;
		${$h_cluster_Mix}{"abundance_in_ref_A"}->{$hap} ||= 0;
		${$h_cluster_Mix}{"z_score_in_ref_B"}->{$hap} ||= 10;
		${$h_cluster_Mix}{"z_score_in_ref_A"}->{$hap} ||= 10;
		${$h_cluster_Mix}{"matched_hap_in_ref_B"}->{$hap} ||= "NA";
		${$h_cluster_Mix}{"matched_hap_in_ref_A"}->{$hap} ||= "NA";
		my $abundance = ${$h_cluster_Mix}{"hapCounts"}->{$hap}/${$h_cluster_Mix}{"total_hap_count"};
		my $abundance_in_A = ${$h_cluster_Mix}{"abundance_in_ref_A"}->{$hap};
		my $abundance_in_B = ${$h_cluster_Mix}{"abundance_in_ref_B"}->{$hap};
		#print "RefA: ", ${$h_cluster_Mix}{"z_score_in_ref_A"}->{$hap}, " RefB:", ${$h_cluster_Mix}{"z_score_in_ref_B"}->{$hap}, "\n";
		#print "RefA: ", ${$h_cluster_Mix}{"abundance_in_ref_A"}->{$hap}, " RefB:", ${$h_cluster_Mix}{"abundance_in_ref_B"}->{$hap}, "\n";
		if(${$h_cluster_Mix}{"z_score_in_ref_A"}->{$hap} - ${$h_cluster_Mix}{"z_score_in_ref_B"}->{$hap} > 0.5){
			$abundance_in_A = 0; #If the difference in Z_score is greater than 0.5, do not count.
		}elsif(${$h_cluster_Mix}{"z_score_in_ref_B"}->{$hap} - ${$h_cluster_Mix}{"z_score_in_ref_A"}->{$hap} > 0.5){
			$abundance_in_B = 0; #If the difference in Z_score is greater than 0.5, do not count.
		}
		#print "$hap\t$abundance, RefA: $abundance_in_A RefB: $abundance_in_B\n";

		$reg_hap->include($abundance, [1, $abundance_in_A, $abundance_in_B]);
		if(($abundance_in_A+$abundance_in_B) > 0.5){ #evidence should be at least 0.5? KZ had 0.001
			$haps_assigned_to_A+=${$h_cluster_Mix}{"hapCounts"}->{$hap}*($abundance_in_A/($abundance_in_A+$abundance_in_B));
			$haps_assigned_to_B+=${$h_cluster_Mix}{"hapCounts"}->{$hap}*($abundance_in_B/($abundance_in_A+$abundance_in_B));
		}else{
			$haps_unassigned+=${$h_cluster_Mix}{"hapCounts"}->{$hap};
		}
	}
	#$reg_hap->print();
	#my $adjrsq = $reg_hap->adjrsq();
	return ($haps_assigned_to_A, $haps_assigned_to_B,$haps_unassigned);
}

#find haplotype groups in B that are similar to A and their relative abundance in A.
sub find_similar_haplotype(){
	my $h_dominantHapInfo_A = shift;
	my $h_dominantHapInfo_B = shift;
	my @similar_haps;
	my @unmatched_haps;
	my @diff_haps;	
	my @haps_A = keys(%{${$h_dominantHapInfo_A}{"haps"}});
	my $expected_distance = $mean_error_rate*length($haps_A[0]);
	foreach my $hap_B (keys(%{${$h_dominantHapInfo_B}{"haps"}})){
		if(${$h_dominantHapInfo_A}{"haps"}->{$hap_B}){
			my $abundance_in_A = sprintf("%6.5f",${$h_dominantHapInfo_A}{"hapCounts"}->{$hap_B}/${$h_dominantHapInfo_A}{"total_hap_count"});
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
		my $abundance_in_A = sprintf("%6.5f",${$h_dominantHapInfo_A}{"hapCounts"}->{$closest_hap}/${$h_dominantHapInfo_A}{"total_hap_count"});
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
sub one_hap_per_cluster(){
	my $h_haplotype_group=shift;
	my %dominantHapInfo;
	foreach my $hapString (keys(%{$h_haplotype_group})){
		$dominantHapInfo{"hapCounts"}->{$hapString}=${$h_haplotype_group}{$hapString};
		$dominantHapInfo{"haps"}->{$hapString}=$hapString;	
		$dominantHapInfo{"total_hap_count"}+=${$h_haplotype_group}{$hapString};
	}
	return%dominantHapInfo;
}

#find all haplotype clusters
sub find_all_haplotype_cluster(){
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
	}while($remainingHaps > 0 && $totalHaps > 0);
	#print "entropy=",$dominantHapInfo{"entropy"},"\n";
	return(%dominantHapInfo);
}

#find the dominant haplotype cluster
sub find_dominant_haplotype_cluster(){
	my $h_haplotype_group=shift;
	my $totalHaps = 0;	
	my $dominant_hap = "";
	my $dominant_hap_count = 0;
	my $dominant_hap_valid_sites=0;
	my @average_methylation;
	foreach my $hapString (keys(%{$h_haplotype_group})){
		#print "hapString:$hapString\t", ${$h_haplotype_group}{$hapString}, "\n";
		$totalHaps+=${$h_haplotype_group}{$hapString};
		my $hap_valid_sites=0;
		for(my $i=0;$i<length($hapString);$i++){
			$average_methylation[$i]=0 if(!$average_methylation[$i]);
			$average_methylation[$i]+=${$h_haplotype_group}{$hapString} if(substr($hapString,$i,1) =~ /[C1]/);
			$hap_valid_sites++ if(substr($hapString,$i,1) ne "N");
		}
		next if($hap_valid_sites<$dominant_hap_valid_sites);
		next if($hap_valid_sites==$dominant_hap_valid_sites && ${$h_haplotype_group}{$hapString} <=$dominant_hap_count);
		$dominant_hap = $hapString;
		$dominant_hap_count = ${$h_haplotype_group}{$hapString};
		$dominant_hap_valid_sites=$hap_valid_sites;
	}
	for(my $i=0;$i<scalar(@average_methylation);$i++){
		$average_methylation[$i]/=$totalHaps;
		$average_methylation[$i] = $mean_error_rate if($average_methylation[$i]<$mean_error_rate);
		$average_methylation[$i] = 1-$mean_error_rate if($average_methylation[$i]>1-$mean_error_rate);
	}
	my $hapsCount_in_cluster=0;
	my @haps_cluster;	
	my $expected_distance = $mean_error_rate*length($dominant_hap);
	foreach my $hapString (keys(%{$h_haplotype_group})){
		next if(!${$h_haplotype_group}{$hapString});
		my $p_value=1;
		my $overlapping_bases=0;
		for(my $i=0; $i<length($hapString); $i++){
			$overlapping_bases++ if(substr($hapString,$i,1) ne "N" && substr($dominant_hap,$i,1) ne "N" );
			next if(substr($hapString,$i,1) eq "N" || 
					substr($dominant_hap,$i,1) eq "N" ||
					substr($hapString,$i,1) eq substr($dominant_hap,$i,1));
			if(substr($dominant_hap,$i,1) =~ /[C1]/){
				$p_value*=1-$average_methylation[$i]+0.0001;			
			}else{
				$p_value*=$average_methylation[$i]+0.0001;	
			}
		}
		next if($p_value < 0.05 || $overlapping_bases<3); #95% or too few overlapping bases
		push(@haps_cluster,$hapString);
		$hapsCount_in_cluster+=${$h_haplotype_group}{$hapString};
	}
	foreach my $hapString (@haps_cluster){
		delete(${$h_haplotype_group}{$hapString});
	}
	#print "####Dominant_hap:$dominant_hap\t$hapsCount_in_cluster\t",join(",",@haps_cluster),"\n" if(!$totalHaps || !$hapsCount_in_cluster);
	return($dominant_hap, join(",",@haps_cluster), $totalHaps, $hapsCount_in_cluster);	
}

#calculate average methylation per target
sub calc_mean_methylation(){
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

sub load_data_files{
	#1) ID
	#2) hapString
	#3) hapCount
	#4) cpgPositions
	#read ref A haplotypes
	my $hap_number=0;
	open(INFILE, "$ref_hap_file_A")||die("Error in opening file $ref_hap_file_A\n");
	while(my $line = <INFILE>){
		chomp($line);
		my @fields = split(/\t/, $line);
		next if(scalar(@fields)<4);
		my $probeID = $fields[0];		
		my $hapString = $fields[1];
		next if(length($hapString)<$minHapLen);
		my $hapCount = $fields[2];
		my @cpgPositions = split(/,/, $fields[3]);
		$hapString =~ s/C/1/g;$hapString =~ s/T/0/g; $hapString =~ s/[AG]/N/g;
		$hapTable{$probeID}->{"hapCounts"}->{"A"}->{$hapString} += $hapCount;
		@{$hapTable{$probeID}->{"cpgPositions"}} = @cpgPositions;
		$hapTable{$probeID}->{"totalHap"}->{"A"}+=$hapCount;
		$hap_number+=$hapCount;
	}
	close(INFILE);
	#print "$hap_number haplotypes loaded for $ref_hap_file_A\n";
	
	#read ref B haplotypes
	$hap_number=0;
	open(INFILE, "$ref_hap_file_B")||die("Error in opening file $ref_hap_file_B\n");
	while(my $line = <INFILE>){
		chomp($line);
		my @fields = split(/\t/, $line);
		next if(scalar(@fields)<4);
		my $probeID = $fields[0];
		my $hapString = $fields[1];
		next if(length($hapString)<3);
		my $hapCount = $fields[2];
		my @cpgPositions = split(/,/, $fields[3]);
		$hapString =~ s/C/1/g;$hapString =~ s/T/0/g;$hapString =~ s/[AG]/N/g;
		$hapTable{$probeID}->{"hapCounts"}->{"B"}->{$hapString} += $hapCount;
		@{$hapTable{$probeID}->{"cpgPositions"}} = @cpgPositions;
		$hapTable{$probeID}->{"totalHap"}->{"B"}+=$hapCount;
		$hap_number+=$hapCount;
	}
	close(INFILE);
	#print "$hap_number haplotypes loaded for $ref_hap_file_B\n";

	#read mix sample haplotypes
	$hap_number=0;
	open(INFILE, "$mix_hap_file")||die("Error in opening file $mix_hap_file\n");
	while(my $line = <INFILE>){
		chomp($line);
		my @fields = split(/\t/, $line);
		next if(scalar(@fields)<4);
		my $probeID = $fields[0];
		my $hapString = $fields[1];
		next if( length($hapString)<3);
		my $hapCount = $fields[2];
		my @cpgPositions = split(/,/, $fields[3]);
		$hapString =~ s/C/1/g;$hapString =~ s/T/0/g;$hapString =~ s/[AG]/N/g;
		$hapTable{$probeID}->{"hapCounts"}->{"mix"}->{$hapString} += $hapCount;		
		@{$hapTable{$probeID}->{"cpgPositions"}} = @cpgPositions;
		$hapTable{$probeID}->{"totalHap"}->{"mix"}+=$hapCount;
		$hap_number+=$hapCount;
	}
	close(INFILE);
	#print "$hap_number haplotypes loaded for $mix_hap_file\n";
}

main(); 
