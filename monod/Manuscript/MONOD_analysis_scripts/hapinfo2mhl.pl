#!/usr/bin/perl -w
# Hapinfo to methylation haplotype load (MHL)
# Run the script to the Hapinfo directory
# Contact: Kun Zhang
# Version 1.3
# Update: 2016-02-29

use strict;
use Cwd;
die &USAGE if @ARGV <1;
my %mch_load_matrix;
my %probe_HMH_samples;
my %hap_count_matrix;
my $hapinfList=shift @ARGV;
open FF,$hapinfList;
chomp(my @hapInfo_files=<FF>);
close FF;

my @sample_list;
foreach my $hapInfo_file(sort @hapInfo_files){
        my @line=split /\//,$hapInfo_file;
	my $sample_name = $line[$#line];
	$sample_name =~ s/.hapInfo.txt//;
	push(@sample_list, $sample_name);
	open(INFILE, "$hapInfo_file") || die("Error in opening $hapInfo_file!");
	while(my $line = <INFILE>){
		chop($line);
		my @fields = split(/\t/, $line);
		next if(scalar(@fields)<4);
		my $probeID = $fields[0];
		my $hapString = $fields[1];
		next if(length($hapString)<1);		
		$hap_count_matrix{$probeID}->{$sample_name}->{$hapString}+=$fields[2];
	}
	close(INFILE);
}

my @unmethylated_haps= ("T", "TT", "TTT", "TTTT", "TTTTT","TTTTTT","TTTTTTT","TTTTTTTT","TTTTTTTTT");
my @methylated_haps  = ("C", "CC", "CCC", "CCCC", "CCCCC","CCCCCC","CCCCCCC","CCCCCCCC","CCCCCCCCC");

foreach my $probeID (keys(%hap_count_matrix)){
	foreach my $sample_name (keys(%{$hap_count_matrix{$probeID}})){
		my %k_mer_counts;
		my $mc_hap_load=0;
		
		foreach my $hapString (keys(%{$hap_count_matrix{$probeID}->{$sample_name}})){
			for(my $word_size = 1; $word_size<=length($hapString); $word_size++){
				next if($word_size>9);
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


print "Probe_id\t", join("\t", sort @sample_list), "\n";
foreach my $probeID (sort keys(%mch_load_matrix)){
	print "$probeID";
	foreach my $sample_name(sort @sample_list){
		$mch_load_matrix{$probeID}->{$sample_name}="NA" if(! defined($mch_load_matrix{$probeID}->{$sample_name}));
		print "\t", $mch_load_matrix{$probeID}->{$sample_name};
	}
	print "\n";
}

sub USAGE{
print "\nperl $0 Hapinfo_File_list > Ouput.txt\n";
print "Just use: ls *hapInfo.txt > Hapinfo_File_list to Get Hapinfo_File_list\n";
}
