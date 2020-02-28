#!/usr/bin/perl -w
use strict;
my $cpg_position_file="/home/kunzhang/HsGenome/hg19/HsGenome19.CpG.positions.txt";
my $samtools = "/home/kunzhang/bin/samtools";
my $probe_info_file= $ARGV[0];
my $fwd_bam_file = $ARGV[1];
my $rev_bam_file = $ARGV[2];
my $UMI_len = $ARGV[3];
my $read_length = 80;
my $use_UMI = 1;
my $PE=1;
$UMI_len = 8 if(!$UMI_len);
die("USAGE: bam2hapInfo.pl probe_info_file fwd_bam rev_bam\n") if(scalar(@ARGV)<3);

my %probeTable;
my %hapInfoTable;
my %binnedProbeTable;
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
my ($reads_parsed, $reads_after_filtering, $UMI_reported);
my %chrSizes=("chr1","249250621","chr2","243199373","chr3","198022430","chr4","191154276",
                                "chr5","180915260","chr6","171115067","chr7","159138663",
                                "chr8","146364022","chr9","141213431","chr10","135534747",
                                "chr11","135006516","chr12","133851895","chr13","115169878",
                                "chr14","107349540","chr15","102531392","chr16","90354753",
                                "chr17","81195210","chr18","78077248","chr19","59128983",
                                "chr20","63025520","chr21","48129895", "chr22","51304566",
                                "chrX","155270560","chrY","59373566","chrM","16571");

open(INFILE, "$probe_info_file")||die("Error in opening $probe_info_file\n");
#forward strand
my $bin_size = 10000;
while(my $line = <INFILE>){
	chop($line);
	my ($chr,$target_start,$target_end,$strand) = split(/[\t\:\-]/, $line);
	my $index = int($target_start/$bin_size);
	if($PE){
		my $probe_id = "$chr:$target_start-$target_end:$strand";
		push(@{$binnedProbeTable{$chr}->{$index}},$probe_id);
		push(@{$binnedProbeTable{$chr}->{$index-1}},$probe_id);
		push(@{$binnedProbeTable{$chr}->{$index+1}},$probe_id);
	}else{
		my $read1_end = $target_start+$read_length;
		my $read2_start = $target_end-$read_length;
		my $probe_id_1 = "$chr:$target_start-$read1_end:$strand";
		my $probe_id_2 = "$chr:$read2_start-$target_end:$strand";
		push(@{$binnedProbeTable{$chr}->{$index}},$probe_id_1);
		push(@{$binnedProbeTable{$chr}->{$index-1}},$probe_id_1);
		push(@{$binnedProbeTable{$chr}->{$index+1}},$probe_id_1);
		push(@{$binnedProbeTable{$chr}->{$index}},$probe_id_2);
		push(@{$binnedProbeTable{$chr}->{$index-1}},$probe_id_2);
		push(@{$binnedProbeTable{$chr}->{$index+1}},$probe_id_2);		
	}
}
close(INFILE);

open(INFILE, "$cpg_position_file")||next;
#forward strand
while(my $line = <INFILE>){
	chop($line);
	my ($chr, $pos) = split(/[\t ]+/, $line);
	my $index = int($pos/$bin_size);
	next if(!$binnedProbeTable{$chr}->{$index});	
	foreach my $probe_id (@{$binnedProbeTable{$chr}->{$index}}){
		my ($chr,$target_start,$target_end,$strand) = split(/[:\-]/, $probe_id);
		if($target_start <= $pos && $pos <=$target_end){
			push(@{$probeTable{$strand}->{$probe_id}->{"CpG_positions"}},$pos);			
		}
	}
}
close(INFILE);

#processing reads mapped to the Watson strand.
foreach my $probe_id (sort keys(%{$probeTable{"W"}})){
	my ($chr,$target_start,$target_end,$strand) = split(/[:\-]/, $probe_id);
	my $region_start = $target_start-30;
	my $region_end = $target_end+30;
	my $cmd = "$samtools view -q 20 $fwd_bam_file $chr:$region_start-$region_end";
	my $ans = `$cmd`;
	my @lines = split(/\n/, $ans);
	$reads_parsed += scalar(@lines);
	next if(scalar(@lines)<1);
	$reads_after_filtering += scalar(@lines);
	my @CpG_positions = sort {$a <=> $b} @{$probeTable{$strand}->{$probe_id}->{"CpG_positions"}};
	my %readHapInfo;	

	#Going through the reads one at a time
	foreach my $line (@lines){
		my @fields = split(/\t/, $line);
		my $seq = $fields[9];
		my $read_start = $fields[3];
		my @read_id_fields = split(/\:/, $fields[0]);
		my $UMI = $read_id_fields[1];	
		#removing the last four elements
		pop(@read_id_fields);pop(@read_id_fields);pop(@read_id_fields);
		pop(@read_id_fields)if($PE);
		my $read_id = join(":", @read_id_fields);	
		foreach my $CpG_position(@CpG_positions){
			my $offset = $CpG_position-$read_start+1;
			next if($offset<0 || $offset >= length($seq));	
			if($use_UMI){
				if(!$readHapInfo{$UMI}->{$read_id}->{$CpG_position}){
					$readHapInfo{$UMI}->{$read_id}->{$CpG_position}=substr($seq,$offset,1);
				}else{
					$readHapInfo{$UMI}->{$read_id}->{$CpG_position}=substr($seq,$offset,1) if(substr($seq,$offset,1) !~/[AGN]/);				
				}
			}else{
				if(!$readHapInfo{$read_id}->{$CpG_position}){
					$readHapInfo{$read_id}->{$CpG_position}=substr($seq,$offset,1);
				}else{
					$readHapInfo{$read_id}->{$CpG_position}=substr($seq,$offset,1) if(substr($seq,$offset,1) !~/[AGN]/);				
				}				
			}			
		}		
	}
		
	if($use_UMI){
		foreach my $UMI (keys(%readHapInfo)){
			next if(!$readHapInfo{$UMI});
			my $hap_string_consensus="";
			foreach my $CpG_position(@CpG_positions){
				my %base_counts=("C",0,"T",0,"N",0);
				foreach my $read_id (keys(%{$readHapInfo{$UMI}})){
					$readHapInfo{$UMI}->{$read_id}->{$CpG_position}="N" if(!$readHapInfo{$UMI}->{$read_id}->{$CpG_position});
					$base_counts{$readHapInfo{$UMI}->{$read_id}->{$CpG_position}}++;
				}
				my $consensus_allele;
				if($base_counts{"C"}>$base_counts{"T"}){
					$consensus_allele = "C";
				}elsif($base_counts{"C"}<$base_counts{"T"}){
					$consensus_allele = "T";
				}else{
					$consensus_allele = "N";
				}
				$hap_string_consensus=$hap_string_consensus. $consensus_allele;
			}
			$hapInfoTable{$probe_id}->{"hap_counts"}->{$hap_string_consensus}++ if($hap_string_consensus);
			$UMI_reported++;
		}
	}else{
		my @read_list = keys(%readHapInfo);	
		for(my $i=0; $i<scalar(@read_list); $i++){		
			$UMI_reported++;
			my $read_id = $read_list[$i];	
			my $hap_string="";
			foreach my $CpG_position(@CpG_positions){
					my $base = $readHapInfo{$read_id}->{$CpG_position} ? $readHapInfo{$read_id}->{$CpG_position} : "N";
					$hap_string=$hap_string.$base;
			}		
			#print "$cmd\t$probe_id\t$read_id\n", join("-",@CpG_positions), "\t", $hap_string, "\n" if($hap_string);
			$hapInfoTable{$probe_id}->{"hap_counts"}->{$hap_string}++ if($hap_string);
		}
	}
	foreach my $hap_string(keys %{$hapInfoTable{$probe_id}->{"hap_counts"}}){
		print $probe_id, "\t", $hap_string, "\t", $hapInfoTable{$probe_id}->{"hap_counts"}->{$hap_string},"\t",join(",",@CpG_positions), "\n";
	}
}

#processing reads mapped to the Crick strand.
foreach my $probe_id (sort keys(%{$probeTable{"C"}})){
	my ($chr,$target_start,$target_end,$strand) = split(/[:\-]/, $probe_id);
	my $region_start_rev = $chrSizes{$chr}-$target_end-30;
	my $region_end_rev = $chrSizes{$chr}-$target_start+30;
	my $cmd = "$samtools view -q 20 $rev_bam_file $chr:$region_start_rev-$region_end_rev";
	my $ans = `$cmd`;
	my @lines = split(/\n/, $ans);
	$reads_parsed += scalar(@lines);
	next if(scalar(@lines)<1);
	$reads_after_filtering += scalar(@lines);
	my @CpG_positions = sort {$a <=> $b} @{$probeTable{$strand}->{$probe_id}->{"CpG_positions"}};
	my %readHapInfo;
	#Going through the reads one at a time
	foreach my $line (@lines){
		my @fields = split(/\t/, $line);
		my $seq = $fields[9];
		my $read_start_rev = $fields[3];
		my @read_id_fields = split(/\:/, $fields[0]);
		#removing the last four elements
		my $UMI = $read_id_fields[1];	
		pop(@read_id_fields);pop(@read_id_fields);pop(@read_id_fields);
		pop(@read_id_fields)if($PE);
		my $read_id = join(":", @read_id_fields);		
		foreach my $CpG_position(@CpG_positions){
			my $offset = $chrSizes{$chr}-$CpG_position-$read_start_rev-1;			
			next if($offset<0 || $offset >= length($seq));				
			#print "Cpg=$CpG_position; offset=$offset; base=",substr($seq,$offset,1),"\t";
			if($use_UMI){
				if(!$readHapInfo{$UMI}->{$read_id}->{$CpG_position}){
					$readHapInfo{$UMI}->{$read_id}->{$CpG_position}=substr($seq,$offset,1);
				}else{
					$readHapInfo{$UMI}->{$read_id}->{$CpG_position}=substr($seq,$offset,1) if(substr($seq,$offset,1) !~/[AGN]/);
				}
			}else{
				if(!$readHapInfo{$read_id}->{$CpG_position}){
					$readHapInfo{$read_id}->{$CpG_position}=substr($seq,$offset,1);
				}else{
					$readHapInfo{$read_id}->{$CpG_position}=substr($seq,$offset,1) if(substr($seq,$offset,1) !~/[AGN]/);
				}
			}
		}	
	}		
		
	if($use_UMI){
		foreach my $UMI (keys(%readHapInfo)){
			next if(!$readHapInfo{$UMI});
			my $hap_string_consensus="";
			foreach my $CpG_position(@CpG_positions){
				my %base_counts=("C",0,"T",0,"N",0);
				foreach my $read_id (keys(%{$readHapInfo{$UMI}})){
					$readHapInfo{$UMI}->{$read_id}->{$CpG_position}="N" if(!$readHapInfo{$UMI}->{$read_id}->{$CpG_position});
					$base_counts{$readHapInfo{$UMI}->{$read_id}->{$CpG_position}}++;
				}
				my $consensus_allele;
				if($base_counts{"C"}>$base_counts{"T"}){
					$consensus_allele = "C";
				}elsif($base_counts{"C"}<$base_counts{"T"}){
					$consensus_allele = "T";
				}else{
					$consensus_allele = "N";
				}
				$hap_string_consensus=$hap_string_consensus. $consensus_allele;
			}
			$hapInfoTable{$probe_id}->{"hap_counts"}->{$hap_string_consensus}++ if($hap_string_consensus);
			$UMI_reported++;	
		}
	}else{
		my @read_list = keys(%readHapInfo);	
		for(my $i=0; $i<scalar(@read_list); $i++){		
			my $read_id = $read_list[$i];	
			my $hap_string="";
			foreach my $CpG_position(@CpG_positions){
					my $base = $readHapInfo{$read_id}->{$CpG_position} ? $readHapInfo{$read_id}->{$CpG_position} : "N";
					$hap_string=$hap_string.$base;
			}		
			#print "$cmd\t$probe_id\t$read_id\n", join("-",@CpG_positions), "\t", $hap_string, "\n" if($hap_string);
			$hapInfoTable{$probe_id}->{"hap_counts"}->{$hap_string}++ if($hap_string);
		}
	}
	
	foreach my $hap_string(keys %{$hapInfoTable{$probe_id}->{"hap_counts"}}){
		print $probe_id, "\t", $hap_string, "\t", $hapInfoTable{$probe_id}->{"hap_counts"}->{$hap_string},"\t",join(",",@CpG_positions), "\n";
	}
}
print STDERR "#reads parsed=$reads_parsed, # reads post-filtering=$reads_after_filtering, #UMI reported=$UMI_reported\n";

# randomly permutate @array in place
sub fisher_yates_shuffle
{
    my $array = shift;
    my $i = @$array;
	return if($i<2);
    while ( --$i )
    {
        my $j = int rand( $i+1 );
        @$array[$i,$j] = @$array[$j,$i];
    }
}

