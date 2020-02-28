#!/usr/bin/perl -w
use strict;
my %mch_load_matrix;
my %probe_HMH_samples;
my %hap_count_matrix;
my $ans = `ls *.hapInfo.txt`;
my @hapInfo_files= split(/[\t\n\r]/, $ans);
my @sample_list;
foreach my $hapInfo_file(@hapInfo_files){
	my $sample_name = $hapInfo_file;
	$sample_name =~ s/.hapInfo.txt//;
	push(@sample_list, $sample_name);
	open(INFILE, "$hapInfo_file") || die("Error in opening $hapInfo_file!");
	while(my $line = <INFILE>){
		chop($line);
		my @fields = split(/\t/, $line);
		next if(scalar(@fields)<4);
		my $probeID = $fields[0];
		my $hapString = $fields[1];
		next if(length($hapString)<4);		
		$hap_count_matrix{$probeID}->{$sample_name}->{$hapString}=$fields[2];
	}
	close(INFILE);
}

my @unmethylated_haps= ("T", "TT", "TTT", "TTTT", "TTTTT");
my @methylated_haps= ("C", "CC", "CCC", "CCCC", "CCCCC");

foreach my $probeID (keys(%hap_count_matrix)){
	foreach my $sample_name (keys(%{$hap_count_matrix{$probeID}})){
		my %k_mer_counts;
		my $mc_hap_load=0;
		
		foreach my $hapString (keys(%{$hap_count_matrix{$probeID}->{$sample_name}})){
			for(my $word_size = 1; $word_size<=length($hapString); $word_size++){
				next if($word_size>5);
				for(my $i=0; $i<=length($hapString)-$word_size; $i++){
					my $sub_hapString = substr($hapString,$i,$word_size);
					next if($sub_hapString =~ /[NAG]/i);
					$k_mer_counts{$word_size}->{$sub_hapString}+=$hap_count_matrix{$probeID}->{$sample_name}->{$hapString};					
				}
			}
		}
		my $norm_factor=0;
		foreach my $word_size (keys(%k_mer_counts)){
			$k_mer_counts{$word_size}->{$unmethylated_haps[$word_size-1]}=0 if(!$k_mer_counts{$word_size}->{$unmethylated_haps[$word_size-1]});
			$k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]}=0 if(!$k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]});
			my $total_count=0;
			foreach my $allele (keys(%{$k_mer_counts{$word_size}})){
				$total_count+=$k_mer_counts{$word_size}->{$allele};
				
			}
			#my $total_count = $k_mer_counts{$word_size}->{$unmethylated_haps[$word_size-1]} + $k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]};
			next if($total_count<1);
			my $mh_fraction = $k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]}/$total_count;
			my $weight = $word_size;
			$mc_hap_load += $weight*$mh_fraction;
			$norm_factor+=$weight;
		}
		next if(!$norm_factor);
		$mc_hap_load/=$norm_factor;
		$mch_load_matrix{$probeID}->{$sample_name}=$mc_hap_load;
	}
}


print "Probe_id\t", join("\t", @sample_list), "\n";
foreach my $probeID (keys(%mch_load_matrix)){
	print "$probeID";
	foreach my $sample_name(@sample_list){
		$mch_load_matrix{$probeID}->{$sample_name}="NA" if(!$mch_load_matrix{$probeID}->{$sample_name});
		print "\t", $mch_load_matrix{$probeID}->{$sample_name};
	}
	print "\n";
}
