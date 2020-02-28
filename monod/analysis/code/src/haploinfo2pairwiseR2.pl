#!/usr/bin/env perl
#
#    Copyright (C) 2014, 2015 Institute for Genomic Medicine (IGM).
#    Portions copyright (C) 2015 Unversity of California, San Diego.
#
#    Author: Kun Zhang <k4zhang@ucsd.edu> and Shicheng Guo <scguo@ucsd.edu>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.


#!/usr/bin/perl -w
use strict;
use Cwd;
use Sort::Array qw/Sort_Table/;

chdir getcwd;

my %hapTable; # probeID=>hapCounts=>sampleID
              #        =>CpgPositions =>$pos
			  #        =>hapCounts
			  # 	   =>totalHap
			  #        =>pairwise_gapb => $gap => $R2

my %LD_matrix;			  
sub main(){
#	my $hapInfo_file = $ARGV[0];
	my $hapInfo_file ="C:\\Users\\shicheng\\Downloads\\STL011LI-01.all_chrs.hapInfo.txt.1000";
	load_hapInfo_files($hapInfo_file);	
	my @all_probe_IDs = sort keys(%hapTable);
	foreach my $probeID (@all_probe_IDs){
		pairwiseR2($probeID);
	}
}			

 
sub load_hapInfo_files(){
	my $hapInfo_file = shift;
	open(INFILE, "$hapInfo_file")||die("Error in opening file $hapInfo_file\n");
	while(my $line = <INFILE>){
		chop($line);
		my @fields = split(/\t/, $line);
		next if(scalar(@fields)<4);
		my $probeID = $fields[0];
		my $hapString = $fields[1];
		next if(length($hapString)<1);
		my $hapCount = $fields[2];
		my @cpgPositions = split(/,/, $fields[3]);
		push(@{$hapTable{$probeID}->{"hapInfo"}->{"Sample"}}, "$hapCount:$hapString:".$fields[3]);
		foreach my $pos(@cpgPositions){
			$hapTable{$probeID}->{"cpgPositionTable"}->{$pos}=1;
		}
	}
	close(INFILE);
	#print "$hap_number haplotypes loaded for $hapInfo_file\n";
	
	my @all_probe_IDs = sort keys(%hapTable);
	foreach my $probeID (@all_probe_IDs){
	    
	    # the following 4 lines can be remove after solve the problem: blank haplotype in the code of bam2haploinfo.pl	
		if(	!$hapTable{$probeID}->{"hapInfo"}->{"Sample"} ){
			delete($hapTable{$probeID});
			next;
		}
		
		my ($chr, $target_start, $target_end) = split(/[:\-]/, $probeID);
		my %cpg_pos_table;
		my %pos2index;
		
		my @unique_pos = sort {$a <=> $b} keys(%{$hapTable{$probeID}->{"cpgPositionTable"}});
		@{$hapTable{$probeID}->{"cpgPositions"}} = @unique_pos;
		my $full_hap_length = scalar(@unique_pos);
		for(my $i=0; $i<$full_hap_length; $i++){
			$pos2index{$unique_pos[$i]}=$i;
		}
		
		foreach my $this_hap_info (@{$hapTable{$probeID}->{"hapInfo"}->{"Sample"}}){
			my ($count, $hapString, $pos_string) = split(/:/, $this_hap_info);
			my @this_hap_positions = split(/,/, $pos_string);
			my @allele_list;
			for(my $i=0; $i<$full_hap_length; $i++){
				push(@allele_list,"N");
			}
			for(my $i=0; $i<scalar(@this_hap_positions); $i++){
				my $pos=$this_hap_positions[$i];
				$allele_list[$pos2index{$pos}]=substr($hapString,$i,1);
				$cpg_pos_table{$pos}+=$count if(substr($hapString,$i,1) =~ /[CT]/);
			}
			my $full_hap_string = join("", @allele_list);
			$hapTable{$probeID}->{"hapCounts"}->{"Sample"}->{$full_hap_string}+=$count;
			$hapTable{$probeID}->{"totalHap"}->{"Sample"}+=$count;
		}		
		
		#print "$probeID\t", join(",", @unique_pos), "\n";
		delete($hapTable{$probeID}->{"hapInfo"});
		delete($hapTable{$probeID}->{"cpgPositionTable"});
		undef(%cpg_pos_table);
		undef(%pos2index);		
	}	
}

sub get_methylation_LD_blocks_greedy(){
	my $probeID = shift;
	my $threshold = shift;
	my $n_loci = scalar(@{$hapTable{$probeID}->{"cpgPositions"}});
	
	my @blocks;
	my %blockSetTable;
	my $start = 0;
	my $i = 0;
	while($i <= $n_loci){
		while(($i < $n_loci) && (&lookupLD($hapTable{$probeID}->{"hapCounts"}->{"Sample"}, $i, $i+1) > $threshold)){
			$i++;
		}
		#print "Search between $start-$i\n";
		push(@blocks, getAllBlockInRegion($hapTable{$probeID}->{"hapCounts"}->{"Sample"}, $start, $i, $threshold));
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
		my $id = sprintf("B%03d", $i+1);
		my @words = split(/:/, $sortedBlocks[$i]);
		$hapTable{$probeID}->{'MLD_blocks'}->{$id}->{'start'} = $words[0];
		$hapTable{$probeID}->{'MLD_blocks'}->{$id}->{'end'} = $words[1];
	}
}



sub pairwiseR2(){
	my $probeID = shift;
	my @cpgpostion=@{$hapTable{$probeID}->{"cpgPositions"}};
	my $n_loci = scalar(@cpgpostion);
	my %hapR2;
	my $i = 0;
	my $r2tmp=0;
	my $num=0;
	while($i < $n_loci-1){
		my $j=$i+1;
		while($j <= $n_loci-1){
		my $gap=$cpgpostion[$j]-$cpgpostion[$i];
		my $R2=lookupLD($hapTable{$probeID}->{"hapCounts"}->{"Sample"}, $i, $j);
		next if ! defined $R2;
		next if $R2 =~/NA/i;
		print "$gap\t$R2\n";
		}
	}
}

sub get_methylation_LD_blocks_Rsquare(){
	my $probeID = shift;
	my $n_loci = scalar(@{$hapTable{$probeID}->{"cpgPositions"}});
	
	my @r2blocks;
	my %blockSetTable;
	my %hapR2;
	my $start = 0;
	my $ivm = 0;
	my $r2tmp=0;
	my $num=0;
	while($ivm < $n_loci-1){
		my $jvm=$ivm+1;
		while($jvm <= $n_loci-1){
		my $r2tmpp=lookupLD($hapTable{$probeID}->{"hapCounts"}->{"Sample"}, $ivm, $jvm);
		$jvm++;
		next if ! defined $r2tmpp;
		next if $r2tmpp =~ m/NA/i;		
		$r2tmp = $r2tmp+ $r2tmpp;
		$num++;
		}
		#print "Search between $start-$i\n";
		$hapR2{$probeID}->{"hapCounts"}->{"Sample"}->{"r2"}=$r2tmp/$num;
		$ivm++;
	}
		my $r2= defined $hapR2{$probeID}->{"hapCounts"}->{"Sample"}->{"r2"} ? $hapR2{$probeID}->{"hapCounts"}->{"Sample"}->{"r2"} : 'NA';
		print "$probeID\t$r2\n";		
}

sub getAllBlockInRegion(){
	my $h_hap_count_table = shift;
	my $start = shift;
	my $end = shift;
	my $threshold = shift;
	my @blocks;
	if($start == $end){
		push(@blocks,$start . ":" . $end);
		return @blocks;
	}
	my ($block_start, $block_end) = findMaxBlockInRegion($h_hap_count_table, $start, $end, $threshold);
	push(@blocks, $block_start. ":" . $block_end);
	if($block_start > $start){
		my @sub_blocks = &getAllBlockInRegion($h_hap_count_table, $start, $block_start-1, $threshold);
		push(@blocks, @sub_blocks);
	}
	if($block_end < $end){
		my @sub_blocks = &getAllBlockInRegion($h_hap_count_table, $block_end+1, $end, $threshold);
		push(@blocks, @sub_blocks);
	}
	return @blocks;
}

sub findMaxBlockInRegion(){
	my $h_hap_count_table = shift;
	my $start = shift;
	my $end = shift;
	my $threshold = shift;
	my $max_block_start=$start;
	my $max_block_end=$start;
	return ($start, $end) if($start == $end);
	for(my $size = $end-$start+1; $size >1; $size--){
		my $good_block=0;
		for(my $i= $start; $i<= $end-$size+1; $i++){
			$good_block = 1;
			for(my $j= $i; $j<$i+$size; $j++){
				for(my $k = $j+1; $k<$i+$size; $k++){
					if((abs($j-$k)<5) && lookupLD($h_hap_count_table, $j, $k) < $threshold) {
						$good_block = 0;
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

sub calculateRsquare(){
	my $h_hap_count_table = shift;
	my $start = shift;
	my $end = shift;
	my ($r2tmp,$num);
		for(my $i= $start; $i<= $end; $i++){
			for(my $j= $i; $j<$end; $j++){
				$r2tmp = lookupLD($h_hap_count_table, $i, $j);
				$num++;
			}
		}
	return($r2tmp/$num);
}

sub lookupLD(){
	my ($h_hap_count_table, $locusA, $locusB) = @_;	
	if (!$LD_matrix{$locusA}->{$locusB}){	
		my %two_locus_hapTable;

		$two_locus_hapTable{"CC"}=0;
		$two_locus_hapTable{"TT"}=0;
		$two_locus_hapTable{"CT"}=0;
		$two_locus_hapTable{"TC"}=0;

		foreach my $hapString (keys(%{$h_hap_count_table})){
			my $two_locus_hap = substr($hapString,$locusA,1).substr($hapString,$locusB,1);
			next if($two_locus_hap =~ /[AGN]/);
			$two_locus_hapTable{$two_locus_hap} += ${$h_hap_count_table}{$hapString};
		}
	
		my ($abs_Dprime,$r2, $abs_d, $abs_Q) = haplotype2LD(\%two_locus_hapTable);
		
		$LD_matrix{$locusA}->{$locusB} = $r2;
		$LD_matrix{$locusB}->{$locusA} = $r2;
	}
	return $LD_matrix{$locusA}->{$locusB};
}


sub haplotype2LD(){
	my $h_two_locus_hapTable = shift;
	my %two_locus_hapTable = %{$h_two_locus_hapTable};
	my @AlleleTable;
	my ($abs_Dprime,$r2, $abs_d, $abs_Q, $mafA, $mafB);
	foreach my $hap (keys(%two_locus_hapTable)){
		if(!$AlleleTable[0][0]){
				$AlleleTable[0][0] = substr($hap,0,1);    # first allele for first locus
		}elsif($AlleleTable[0][0] ne substr($hap,0,1)){   
				$AlleleTable[0][1] = substr($hap,0,1);    # second allele for first locus
		}
		if(!$AlleleTable[1][0]){
				$AlleleTable[1][0] = substr($hap,1,1);    # first allele for second locus
		}elsif($AlleleTable[1][0] ne substr($hap,1,1)){
				$AlleleTable[1][1] = substr($hap,1,1);    # second allele for second locus
		}
	}

	foreach my $alleleA (@{$AlleleTable[0]}){
		foreach my $alleleB (@{$AlleleTable[1]}){
			if(!$two_locus_hapTable{$alleleA.$alleleB}){
				$two_locus_hapTable{$alleleA.$alleleB} = 0;   # initial the counts of missing allele to zero
			}
		}
	}

	my @AlleleFreq;
	$AlleleFreq[0][0] = $two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][0]} + $two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][1]};  # TC+TT or CT+CC: counts of 1st allele in 1st locus
	$AlleleFreq[0][1] = $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][0]} + $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][1]};  # CT+CC or TC+TT: counts of 2st allele in 1st locus
	$AlleleFreq[1][0] = $two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][0]} + $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][0]};  # CT+TT or TC+CC: counts of 1st allele in 2st locus
	$AlleleFreq[1][1] = $two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][1]} + $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][1]};  # TC+CC or CT+TT: counts of 2st allele in 2st locus
#	return ("NA","NA","NA","NA","NA","NA") if($AlleleFreq[0][0] + $AlleleFreq[0][1] ==0 || $AlleleFreq[1][0] + $AlleleFreq[1][1] == 0);           # skip if any loci have no value

	$mafA = $AlleleFreq[0][0] < $AlleleFreq[0][1] ? $AlleleFreq[0][0] : $AlleleFreq[0][1];                                                        # minor allele frequency for 1st locus
	$mafA /= $AlleleFreq[0][0] + $AlleleFreq[0][1]+0.0000001;
	$mafB = $AlleleFreq[1][0] < $AlleleFreq[1][1] ? $AlleleFreq[1][0] : $AlleleFreq[1][1];                                                        # minor allele frequency for 2st locus
	$mafB /= $AlleleFreq[1][0] + $AlleleFreq[1][1]+0.0000001;
#	return ("NA","1","NA","NA", $mafA, $mafB) if($mafA ==0 && $mafB==0);
	
	# $D=PAB-PA*PB
	# Dprime=$D/$Dmax
	# $Dmax= $D<0 ? min(PAB,(1-PA)*(1-PB) : min(PA*(1-PB),(1-PA)*PB)
	
	my $D = $two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][0]} * $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][1]} 
	      - $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][0]} * $two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][1]};
	my $Dmax = $D > 0 ? ($AlleleFreq[0][0]*$AlleleFreq[1][1] < $AlleleFreq[0][1]*$AlleleFreq[1][0] ? $AlleleFreq[0][0]*$AlleleFreq[1][1] :$AlleleFreq[0][1]*$AlleleFreq[1][0]) : 
	                    ($AlleleFreq[0][0]*$AlleleFreq[1][0] < $AlleleFreq[0][1]*$AlleleFreq[1][1] ? $AlleleFreq[0][0]*$AlleleFreq[1][0] : $AlleleFreq[0][1]*$AlleleFreq[1][1]);
	if($D == 0.0 ){
		$abs_Dprime = ($Dmax == 0.0)? 1.0 :0.0;
		$r2 = ($AlleleFreq[0][0]*$AlleleFreq[0][1]*$AlleleFreq[1][0]*$AlleleFreq[1][1]) == 0.0 ? 1.0 : 0.0;
		$abs_d = $AlleleFreq[1][0]*$AlleleFreq[1][1] == 0.0 ? 1.0 : 0.0;
		$abs_Q = ($two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][0]} * $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][1]}
				+ $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][0]} * $two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][1]} == 0.0) ? 1.0 : 0.0;
	}else{
		$abs_Dprime = abs($D/$Dmax);
		$r2 = $D*$D/($AlleleFreq[0][0]*$AlleleFreq[0][1]*$AlleleFreq[1][0]*$AlleleFreq[1][1]);
		$abs_d = abs($D/($AlleleFreq[1][0]*$AlleleFreq[1][1]));
		$abs_Q = abs($D/($two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][0]} * $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][1]}
		     + $two_locus_hapTable{$AlleleTable[0][1].$AlleleTable[1][0]} * $two_locus_hapTable{$AlleleTable[0][0].$AlleleTable[1][1]}));
	}
	return ($abs_Dprime,$r2, $abs_d, $abs_Q, $mafA, $mafB);
}

main();