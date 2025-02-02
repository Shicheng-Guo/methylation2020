#!/usr/bin/perl -w
use strict;
my $cpg_position_file="/home/kunzhang/HsGenome/hg19/HsGenome19.CpG.positions.txt";
my $samtools = "/home/kunzhang/bin/samtools";
my $target_list_file= $ARGV[0];
my $merged_bam_file = $ARGV[1];
die("USAGE: bam2hapInfo.pl target_list_file merged_bam\n") if(scalar(@ARGV)<2);

my $target_flanking_len = 80;
my %targetTable;
my %hapInfoTable;
my %binnedtargetTable;
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
my ($reads_parsed, $reads_after_filtering);
my %chrSizes=("chr1","249250621","chr2","243199373","chr3","198022430","chr4","191154276",
                                "chr5","180915260","chr6","171115067","chr7","159138663",
                                "chr8","146364022","chr9","141213431","chr10","135534747",
                                "chr11","135006516","chr12","133851895","chr13","115169878",
                                "chr14","107349540","chr15","102531392","chr16","90354753",
                                "chr17","81195210","chr18","78077248","chr19","59128983",
                                "chr20","63025520","chr21","48129895", "chr22","51304566",
                                "chrX","155270560","chrY","59373566","chrM","16571");

open(INFILE, "$target_list_file")||die("Error in opening $target_list_file\n");
#forward strand
my $bin_size = 10000;
while(my $line = <INFILE>){
	chop($line);
	my ($chr,$target_start,$target_end,$id) = split(/[\t ]+/, $line);
	my $index = int($target_start/$bin_size);
	my $target_id = "$chr:$target_start-$target_end";
	push(@{$binnedtargetTable{$chr}->{$index}},$target_id);
	push(@{$binnedtargetTable{$chr}->{$index-1}},$target_id);
	push(@{$binnedtargetTable{$chr}->{$index+1}},$target_id);
}
close(INFILE);

open(INFILE, "$cpg_position_file")||next;
while(my $line = <INFILE>){
	chop($line);
	my ($chr, $pos) = split(/[\t ]+/, $line);
	my $index = int($pos/$bin_size);
	next if(!$binnedtargetTable{$chr}->{$index});	
	foreach my $target_id (@{$binnedtargetTable{$chr}->{$index}}){
		my ($chr,$target_start,$target_end) = split(/[:\-]/, $target_id);
		if($target_start <= $pos && $pos <=$target_end){
			push(@{$targetTable{$target_id}->{"CpG_positions"}},$pos);			
		}
	}
}
close(INFILE);

foreach my $target_id (sort keys(%targetTable)){
	my ($chr,$target_start,$target_end) = split(/[:\-]/, $target_id);
	my $region_start = $target_start-$target_flanking_len;
	my $region_end = $target_end+$target_flanking_len;
	my $cmd = "$samtools view -q 20 $merged_bam_file $chr:$region_start-$region_end";
	my $ans = `$cmd`;
	my @lines = split(/\n/, $ans);
	#print "##$cmd\n$ans\n\n";
	next if(scalar(@lines)<1 || !$targetTable{$target_id}->{"CpG_positions"});
	$reads_parsed += scalar(@lines);
	$reads_after_filtering += scalar(@lines);
	my @CpG_positions = sort {$a <=> $b} @{$targetTable{$target_id}->{"CpG_positions"}};
	my %readHapInfo;	

	#print "##$cmd\n";
	#Going through the reads one at a time
	foreach my $line (@lines){
		my @fields = split(/\t/, $line);
		my $read_strand = $fields[1];		
		my $seq = $fields[9];
		my $qual_string = $fields[10];
		my $read_start = $fields[3];
		my $CIGAR = $fields[5];
		my $read_length = length($seq);
		next if($CIGAR =~ /[ID]/);
		if($CIGAR=~ /S/){
			my ($clip_len, @others) = split(/S/, $CIGAR);
			next if(length($clip_len)>2);
			$seq=substr($seq,$clip_len,$read_length-$clip_len);		
			$read_length -= $clip_len;			
		}elsif($seq =~ /^CG/){
			$seq=substr($seq,2,$read_length-2);		
			$read_start += 2;
			$read_length -= 2;
		}
		my @read_fields = split(/[\:\#\_]+/, $fields[0]);
		my $read_id = scalar(@read_fields)>8 ? join("-", @read_fields[0..6]) : join("-", @read_fields[0..4]);	
		$read_fields[7] = "NA" if(!$read_fields[7]);
		my $UMI_tag = scalar(@read_fields)>8 ? $read_fields[11] : $read_fields[7];
		$UMI_tag = $UMI_tag . ":".$fields[3];
		#print "##$read_id\t$UMI_tag\t$CIGAR\t$read_start\t$read_strand\t$seq\n";
		
		next if(!$UMI_tag);
		foreach my $CpG_position(@CpG_positions){
			my $offset = $CpG_position-$read_start+1;
			next if($offset<0 || $offset >= length($seq));	
			if($read_strand<16){
				next if($readHapInfo{$UMI_tag}->{$CpG_position}->{"qual"} && 
					$readHapInfo{$UMI_tag}->{$CpG_position}->{"qual"}>=ord(substr($qual_string,$offset,1)));
				$readHapInfo{$UMI_tag}->{$CpG_position}->{"base"}=substr($seq,$offset,1);
				$readHapInfo{$UMI_tag}->{$CpG_position}->{"qual"}=ord(substr($qual_string,$offset,1));
			}else{
				next if($readHapInfo{$UMI_tag}->{$CpG_position}->{"qual"} &&
					$readHapInfo{$UMI_tag}->{$CpG_position}->{"qual"}>=ord(substr($qual_string,$offset+1,1)));
				$readHapInfo{$UMI_tag}->{$CpG_position}->{"base"}=$rcTable{substr($seq,$offset+1,1)};	
				$readHapInfo{$UMI_tag}->{$CpG_position}->{"qual"}=ord(substr($qual_string,$offset+1,1));						
			}	
			#print "$offset\t", $readHapInfo{$read_id}->{$CpG_position}, "\n";
		}	
	}		
	my @UMI_list = keys(%readHapInfo);	
	for(my $i=0; $i<scalar(@UMI_list); $i++){		
			my $UMI_tag = $UMI_list[$i];	
			my $hap_string="";
			foreach my $CpG_position(@CpG_positions){
					my $base = $readHapInfo{$UMI_tag}->{$CpG_position}->{"base"} ? $readHapInfo{$UMI_tag}->{$CpG_position}->{"base"} : "N";
					$hap_string=$hap_string.$base;
			}		
			#print "$target_id\t$UMI_tag\n", join("-",@CpG_positions), "\t", $hap_string, "\n" if($hap_string);
			$hapInfoTable{$target_id}->{"hap_counts"}->{$hap_string}++ if($hap_string);
	}	
	foreach my $hap_string(keys %{$hapInfoTable{$target_id}->{"hap_counts"}}){
		my ($valid_start, $valid_end);
		my @valid_positions;
		my $valid_hap="";
		for(my $i=0; $i<length($hap_string); $i++){
			my $allele = substr($hap_string,$i,1);
			next if($allele eq "N");
			$valid_hap = $valid_hap.$allele;
			push(@valid_positions,$CpG_positions[$i]);
		}
		next if(scalar(@valid_positions)<3);
		print $target_id,"\t$valid_hap\t", $hapInfoTable{$target_id}->{"hap_counts"}->{$hap_string},"\t",join(",",@valid_positions), "\n";
	}
}
	
