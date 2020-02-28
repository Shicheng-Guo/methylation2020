#!/usr/bin/perl -w

use strict;
use warnings;
use Cwd;
use Getopt::Long;
use Pod::Usage;
chdir getcwd;

my ($hapinfList,$bedregion,$output,$help,$version)=&process_command_line();

my %mch_load_matrix;
my %probe_HMH_samples;
my %hap_count_matrix;

open FF,$hapinfList || die "cannot open $hapinfList\n";
chomp(my @hapInfo_files=<FF>);
close FF;

open F,"<$bedregion" || die "cannot open $bedregion\n";;
my %hash;
my $window=5000;
while(<F>){
        chomp;
        next if /^\s+$/;
        my ($chr,$start,$end,$id,$gene,$block)=split/\s+/;
        my $bin1=int($start/$window);
        my $bin2=int($end/$window);
        foreach my $i($bin1..$bin2){
        push @{$hash{$chr}{$i-1}},"$chr:$start-$end";
        push @{$hash{$chr}{$i}},"$chr:$start-$end";
        push @{$hash{$chr}{$i+1}},"$chr:$start-$end";
	}
}
close F;

my @sample_list;

foreach my $hapInfo_file(@hapInfo_files){
    my ($sample_name,undef)=reverse(split /[\/|\t|\s+]/,$hapInfo_file);
	# print "$sample_name\n";  # debug
	push(@sample_list, $sample_name);
	open(INFILE, "$sample_name") || die("Error in opening $sample_name!");
	while(my $line = <INFILE>){
		chomp($line);
		my @fields = split(/\s+/, $line);
		next if(scalar(@fields)<4);
		my $hapString = $fields[1];
		next if(length($hapString)<1);	
		my $probeID = $fields[0];
		my ($chr,undef,undef)=split/[:-]/,$probeID;
		my @tmp=split/,/,$fields[3];
		my $bin=int($tmp[0]/$window);
		if(defined $hash{$chr}{$bin}){
			my $readstart=$tmp[0];
			my $readend=$tmp[$#tmp];
			# print "$line\t";         # debug
			foreach my $bed(sort @{$hash{$chr}{$bin}}){
				# print "$bed\t";    # debug
				my $HapString;
				my ($chr,$start,$end)=split/[:-]/,$bed;
				if($readstart>=$start){
				my @pos;
				foreach my $i(0..$#tmp){
				push @pos,$i if ($tmp[$i]<=$end);
				}
				if(scalar(@pos)==1){
				$HapString=substr($hapString,$pos[0],1);
				$hap_count_matrix{$bed}->{$sample_name}->{$HapString}+=$fields[2];
				}elsif(scalar(@pos)>1){
				$HapString=substr($hapString,$pos[0],scalar(@pos));
				$hap_count_matrix{$bed}->{$sample_name}->{$HapString}+=$fields[2];
				}
				# my $len=scalar(@pos);   # debug
				# print "$len::$HapString\t";  # debug
				}else{
				my @pos;
				foreach my $i(0..$#tmp){
				push @pos,$i if ($tmp[$i]<=$end and $tmp[$i]>=$start);	
				print "$i," if ($tmp[$i]<=$end and $tmp[$i]>=$start);	
				}
				if(scalar(@pos)==1){
				$HapString=substr($hapString,$pos[0],1);
				$hap_count_matrix{$bed}->{$sample_name}->{$HapString}+=$fields[2];
				}elsif(scalar(@pos)>1){
				$HapString=substr($hapString,$pos[0],scalar(@pos));
				$hap_count_matrix{$bed}->{$sample_name}->{$HapString}+=$fields[2];
				}
				# my $len=scalar(@pos);        # debug
				# print "$len::$HapString\t";  # debug
				}		
			}
				# print "\n";                  # debug
		}
	}
	close(INFILE);
}

my @unmethylated_haps= ("T"x1,"T"x2,"T"x3,"T"x4,"T"x5,"T"x6,"T"x7,"T"x8,"T"x9,"T"x10,"T"x11,"T"x12,"T"x13,"T"x14);
my @methylated_haps  = ("C"x1,"C"x2,"C"x3,"C"x4,"C"x5,"C"x6,"C"x7,"C"x8,"C"x9,"C"x10,"C"x11,"C"x12,"C"x13,"C"x14);

# metyalation hapltoype and methylations 
# rnaseq while nameside %% jim
# regulationsd
open OUT,">$output.mhl.txt";
foreach my $probeID (sort keys(%hap_count_matrix)){
	# print "$probeID\n";          # debug
	foreach my $sample_name (sort keys(%{$hap_count_matrix{$probeID}})){
		# print "$sample_name\t";  # debug
		my %k_mer_counts;
		my $mc_hap_load=0;
		foreach my $hapString (sort keys(%{$hap_count_matrix{$probeID}->{$sample_name}})){
			for(my $word_size = 1; $word_size<=length($hapString); $word_size++){
				next if($word_size>10);
				for(my $i=0; $i<=length($hapString)-$word_size; $i++){
					my $sub_hapString = substr($hapString,$i,$word_size);
					next if($sub_hapString =~ /[NAG]/i);
					$k_mer_counts{$word_size}->{$sub_hapString}+=$hap_count_matrix{$probeID}->{$sample_name}->{$hapString};					
				}
			}
		}
		my $norm_factor=0;
		foreach my $word_size (sort keys(%k_mer_counts)){
			$k_mer_counts{$word_size}->{$unmethylated_haps[$word_size-1]}=0 if(!$k_mer_counts{$word_size}->{$unmethylated_haps[$word_size-1]});
			$k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]}=0 if(!$k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]});
			my $total_count=0;
			foreach my $allele (sort keys(%{$k_mer_counts{$word_size}})){
				$total_count+=$k_mer_counts{$word_size}->{$allele};
			}
			next if($total_count<1);
			my $mh_fraction = $k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]}/$total_count;
			my $weight = $word_size;
			$mc_hap_load += $weight*$mh_fraction;
			$norm_factor+=$weight;
			# print "($weight:$mh_fraction:$mc_hap_load)\t";     # debug
			
		}
		    # print "\n";                                        # debug
		next if(!$norm_factor);
		$mc_hap_load/=$norm_factor;
		$mch_load_matrix{$probeID}->{$sample_name}=$mc_hap_load;
	}
}

print OUT "Probe_id\t", join("\t", sort @sample_list), "\n";
foreach my $probeID (sort keys(%mch_load_matrix)){
	print OUT "$probeID";
	foreach my $sample_name(sort @sample_list){
		$mch_load_matrix{$probeID}->{$sample_name}="NA" if(! defined($mch_load_matrix{$probeID}->{$sample_name}));
		print OUT "\t", $mch_load_matrix{$probeID}->{$sample_name};
        # my $mhl=sprintf("%.3f",$mch_load_matrix{$probeID}->{$sample_name});
        # print "\t",$mhl;
	}
	print OUT "\n";
}

sub process_command_line{
	my $hapinfList;
	my $bedregion;
	my $output;
	my $help;
	my $version;	
    my $command_line= GetOptions ("haplist=s"=> \$hapinfList,
              			     "bed=s" => \$bedregion,
	          		     "output=s" => \$output,
	          		     "help" => \$help,
	          		      "version" => \$version,
	          			);
	if ($help){
    	print_helpfile();
    	exit;
    	}  
    	if ($version){
    	version();
    	exit;
    	}  
    	unless(defined $hapinfList and defined $bedregion){
    	print_helpfile();
    	die "\n--haplist, --bed and --output should be provided!\n";
    	}
	 return($hapinfList,$bedregion,$output,$help,$version);
}
 
sub print_helpfile{

	print<<"HOW_TO";

DESCRIPTION

USAGE: hapinfo2mhl --haplist hapinfolist --bed interest.bed --output project.mhl.txt 

Last edited on 05 May 2017 Shicheng Guo <Shihcheng.Guo\@gmail.com>.
	
ARGUMENTS:
		
 --haplist          File list of the haplotype files (hapInfo.txt). All or parts of hapinfo file could be collected. 
	
 --bed              Genomic regions of interesting, such as MHB, Enhancer or your custom defined regions.

 --output           Prefix of the output files. It will be add mhl.txt automatically. 
	
HOW_TO
exit;
	
}

sub version{
	print<<"VERSION";
	
          bam2hapinfo - Smart tools to calculate methylation haplotypes load from methylation haplotype files(hapinfo).

                               bam2hapinfo Version: "bam2hapinfo version: 1.02"
                           
                 Copyright 2010-15 Kun Zhang <kzhang\@eng.ucsd.edu >, University of California, San Diego

			                      Software Maintainer: <shg047\@ucsd.edu>

VERSION
exit;	
}



