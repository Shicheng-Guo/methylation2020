#!/usr/bin/perl -w
# this script reads in a haploInfo file and output the rsq value for each pair of cg sites

use strict;
use warnings;
use Sort::Array qw/Sort_Table/;

my %cgTable;
my %RsqPvalTable;
my $min_coverage = 10;
my $min_rsq = 0.3;

my %rcTable;
$rcTable{'A'}='T';
$rcTable{'T'}='A';
$rcTable{'G'}='C';
$rcTable{'C'}='G';
$rcTable{'N'}='N';
$rcTable{'D'}='D';

printUsage() if(!$ARGV[0]);

sub main{
	readCGTable($ARGV[1]);
	my $last_targets = "NA";
	open(RSQ_OUT, ">$ARGV[0]") || die("error writing to file $ARGV[0]\n");
	print RSQ_OUT "index\tRsq\tmm:mu:um:uu\ttotalHaps\n";
	my %RsqTable;
	my ($chr, $start, $end);
	while(my $line = <STDIN>){
		chomp($line);
		#chr10:10000873-10001472	CCC	1	10001056,10001082,10001168
		my @tmp =  split /\t/, $line;
		if($last_targets ne "NA" && $last_targets ne $tmp[0]){
			# get LD blocks for the last targets
			($chr, $start, $end) = split /[:-]/, $last_targets;
			my %posList;
			foreach my $pairs(keys %RsqTable){
				($chr, $start, $end) = split /[:-]/, $pairs;
				my ($rsq, $pval, $counts, $totalHaps) = split "\t", $RsqTable{$pairs};
				next if($totalHaps < $min_coverage);
				$posList{$start} = 1;
				$posList{$end} = 1;
				print RSQ_OUT $pairs, "\t", $RsqTable{$pairs},"\n";
			}
			my @sorted_list = sort {$a<=>$b} keys %posList;
			if(scalar(@sorted_list) >= 3){
				#getDI(\%RsqTable, \@sorted_list, $chr)
				my %blockSetTable = getBlocksGreedy(\%RsqTable, \@sorted_list, $chr, $min_rsq);
				reportBlocks($ARGV[0].".LDblock.txt", \%blockSetTable);
			}
			undef(%posList);
			undef(%RsqTable);
		}
		$last_targets = $tmp[0];
		print $last_targets, "\n";
		($chr, $start, $end) = split /[:-]/, $tmp[0];
		my ($hapString, $hapCount) = ($tmp[1], $tmp[2]);
		my @cgPos = split ",", $tmp[3];
		next if(length($tmp[1]) ne scalar(@cgPos));
		my %indices;
		for(my $i = 0; $i < scalar(@cgPos); $i++){
			my $pos = $cgPos[$i];
			next if(!$cgTable{$chr.":".$pos}); # is pos an acceptable CG site?
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
		for(my $i = 0; $i < scalar(@ordered_pos); $i++){
			for(my $j = $i+1; $j < scalar(@ordered_pos); $j++){
				my $value_key = $chr.":".$ordered_pos[$i]."-".$ordered_pos[$j];
				my $value_call = $indices{$ordered_pos[$i]}.$indices{$ordered_pos[$j]};
				#$value_call =~ tr/T/U/;
				#$value_call =~ tr/C/M/;
				my ($mm, $mu, $um, $uu) = (0,0,0,0);
				my ($rsq, $pval, $counts, $totalHaps) = ("NA", "NA", "NA", 0);
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
					my $mafA = ($um + $uu) < ($mm + $mu) ? $um + $uu : $mm + $mu;
					my $mafB = ($mu + $uu) < ($mm + $um) ? $mu + $uu : $mm + $um;
					if($mafA ne 0 and $mafB ne 0){
						($rsq, $pval) = updateRsq($counts);
						$RsqTable{$value_key} = "$rsq\t$pval\t$counts\t$totalHaps";
					}else{
						$RsqTable{$value_key} = "0\tMAF_0\t$counts\t$totalHaps";
					}
				}else{
					$RsqTable{$value_key} = "$rsq\t$pval\t$counts\t$totalHaps";
				}		
			}
		}
		undef(@ordered_pos);
		undef(%indices);
	}
	($chr, $start, $end) = split /[:-]/, $last_targets;
	my %posList;
	foreach my $pairs(keys %RsqTable){
		($chr, $start, $end) = split /[:-]/, $pairs;
		my ($rsq, $pval, $counts, $totalHaps) = split "\t", $RsqTable{$pairs};
		next if($totalHaps < $min_coverage);
		$posList{$start} = 1;
		$posList{$end} = 1;
		print RSQ_OUT $pairs, "\t", $RsqTable{$pairs},"\n";
	}
	my @sorted_list = sort {$a<=>$b} keys %posList;
		if(scalar(@sorted_list) >= 3){
		#getDI(\%RsqTable, \@sorted_list, $chr)
		my %blockSetTable = getBlocksGreedy(\%RsqTable, \@sorted_list, $chr, $min_rsq);
		reportBlocks($ARGV[0].".LDblock.txt", \%blockSetTable);
	}
	undef(%posList);
	undef(%RsqTable);
	undef(%cgTable);
	close(RSQ_OUT);
}

sub readCGTable{
        my $file = shift;
        open(IN, "$file") || die("Error opening $file\n");
        while(my $line = <IN>){
                chomp($line);
                my @tmp = split /\t/, $line;
                $tmp[0] =~ s/:W//;
                $tmp[1]--;
                $cgTable{$tmp[0].":".$tmp[1]} = 1;
        }
        close(IN);
}

sub reportBlocks{
	my $filename = shift;
	my $h_blocks = shift;
	my %blocktable = %{$h_blocks};
	open(OUTFILE, ">>$filename") || die("Error opening $filename\n");
	foreach my $id (sort keys(%blocktable)) {
		my ($chr, $n) = split ":", $id;
		next if($blocktable{$id}->{'end'} - $blocktable{$id}->{'start'} < 100);
		print OUTFILE $chr, "\t", $blocktable{$id}->{'start'}, "\t", $blocktable{$id}->{'end'}, "\t", $blocktable{$id}->{'numCG'}, "\n";
	}
	close(OUTFILE);
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
	my $rsq = 0;
	$rsq = ($D*$D)/($p1*$p2*$q1*$q2);
	my $out_rsq = sprintf("%4.3f", $rsq);
	my $out_pval;
	my ($mm_s, $mu_s, $um_s, $uu_s) = (int($min_coverage*$mm), int($min_coverage*$mu), int($min_coverage*$um), int($min_coverage*$uu));
	my $counts_s = "$mm_s:$mu_s:$um_s:$uu_s";
	if(defined($RsqPvalTable{$counts_s})){
		return $out_rsq, $RsqPvalTable{$counts_s};
	}
	# begins permutation:
	my $num_perm = 100;
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
	my $RsqTable_ref = shift;
	my $chr = shift;
	my $pos1 = shift;
	my $pos2 = shift;
	my $rsq;
	if(!${$RsqTable_ref}{$chr.":".$pos1."-".$pos2} && $pos2-$pos1 <= 100){
		$rsq = 1;
	}elsif(${$RsqTable_ref}{$chr.":".$pos1."-".$pos2}){
		my ($rsq_out, $pval, $allele_count, $totalHaps) = split "\t", ${$RsqTable_ref}{$chr.":".$pos1."-".$pos2};
		if($pval eq "MAF_0" || $totalHaps < $min_coverage){
			$rsq = 1;
		}else{
			$rsq = $rsq_out;
		}
	}elsif($pos2-$pos1 > 100){
		$rsq = 0;
	}else{
		$rsq = 0;
	}
	return $rsq;
}

sub getDI{
	my $RsqTable_ref = shift;
	my $sorted_list = shift;
	my $chr = shift;
	my $start = 0;
	my $end = scalar(@{$sorted_list})-1;
	my $i = $start;
	my $max_distance = 200000; #2MB
	my %countTableDown;
	my %countTableUp;
	$countTableDown{$i} = "0\t0\tNA";
	while($i < $end){
		my $j = $i+1;
		my $sum = 0;
		my $count = 0;
		while($j < $end && @{$sorted_list}[$j] - @{$sorted_list}[$i] < $max_distance){
			$sum += lookupLD($RsqTable_ref, $chr, @{$sorted_list}[$i], @{$sorted_list}[$j]);
			$count++;
			if(!$countTableDown{$j}){
				$countTableDown{$j} = "$sum\t$count\t$i";
			}
			$j++;
		}
		$countTableUp{$i} = "$sum\t$count\t$j";
		if(!$countTableDown{$j}){ # choose the farthest i from j for downstream
			$countTableDown{$j} = "$sum\t$count\t$i";
		}
		$i++;
	}
	$countTableUp{$i} = "0\t0\tNA";
	# Calculate DI:
	# A = average R^2 upstream 2MB
	# B = average R^2 downstream 2MB
	# E = expected R^2 (A+B)/2
	my $file = $ARGV[0]. ".DI.txt";
	open(OUT, ">$file") || die("Error writing to $file\n");
	print OUT "Position\tA\tB\tE\tDI\n";
	$i = $start;
	my ($c_sum, $c_count, $c_j);
	while($i <= $end){
		my $c_pos = @{$sorted_list}[$i];
		my ($A_ave, $B_ave, $E_val, $DI) = (0,0,0,0);
		if(!$countTableUp{$i} || !$countTableDown{$i}){
			print OUT "Can't calculate DI for position $c_pos, missing upstream and downstream counts!\n";
		
		}else{
			($c_sum, $c_count, $c_j) = split "\t", $countTableUp{$i};
			$A_ave = $c_sum/$c_count if($c_count > 0);
			($c_sum, $c_count, $c_j) = split "\t", $countTableDown{$i};
			$B_ave = $c_sum/$c_count if($c_count > 0);
			$E_val = ($A_ave + $B_ave)/2;
			$DI = ( ($B_ave - $A_ave)/abs($B_ave - $A_ave) )*( ($A_ave - $E_val)**2/$E_val + ($B_ave - $E_val)**2/$E_val ) if($E_val ne 0 || abs($B_ave - $A_ave) ne 0);
			print OUT $chr, ":", $c_pos, "\t", $A_ave, "\t", $B_ave, "\t", $E_val, "\t", $DI, "\n";
		}
		$i++;
	}
	close(OUT);
}

sub getBlocksGreedy{
	my $RsqTable_ref = shift;
	my $sorted_list = shift;
	my $chr = shift;
	my $start = 0;
	my $end = scalar(@{$sorted_list}) - 1;
	my $threshold = shift;
	my $i = $start;
	my @blocks;
	my %blockSetTable;
	#print "Max i is $end for $chr\n";
	while($i <= $end){
		while(($i < $end) and (lookupLD($RsqTable_ref, $chr, @{$sorted_list}[$i], @{$sorted_list}[$i+1]) > $threshold)){
			$i++;
		}
		#print $i, "\n";
		push(@blocks, getAllBlockInRegion($RsqTable_ref, $sorted_list, $chr, $start, $i, $threshold));
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
		my $id = sprintf("$chr:%05d", $i+1);
		my @words = split(/:/, $sortedBlocks[$i]);
		$blockSetTable{$id}->{'start'} = $words[0];
		$blockSetTable{$id}->{'end'} = $words[1];
		$blockSetTable{$id}->{'numCG'} = $words[2];
	}
	return %blockSetTable;
}

sub getAllBlockInRegion{
	my $RsqTable_ref = shift;
	my $sorted_list = shift;
	my $chr = shift;
	my $start = shift;
	my $end = shift;
	my $threshold = shift;
	my @blocks;
	my $numCG;
	if($start == $end){
		$numCG = $end - $start + 1;
		push(@blocks, @{$sorted_list}[$start] . ":" . @{$sorted_list}[$end] . ":" . $numCG);
		return @blocks;
	}
	my ($block_start, $block_end) = findMaxBlockInRegion($RsqTable_ref, $sorted_list, $chr, $start, $end, $threshold);
	$numCG = $block_end - $block_start + 1;
	push(@blocks, @{$sorted_list}[$block_start] . ":" . @{$sorted_list}[$block_end] . ":" . $numCG) if($block_end - $block_start + 1 > 10);
	if($block_start > $start){
		my @sub_blocks = getAllBlockInRegion($RsqTable_ref, $sorted_list, $chr, $start, $block_start-1, $threshold);
		push(@blocks, @sub_blocks);
	}
	if($block_end < $end){
		my @sub_blocks = getAllBlockInRegion($RsqTable_ref, $sorted_list, $chr, $block_end+1, $end, $threshold);
		push(@blocks, @sub_blocks);
	}
	return @blocks;
}

sub findMaxBlockInRegion(){
	my $RsqTable_ref = shift;
	my $sorted_list = shift;
	my $chr = shift;
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
					if(lookupLD($RsqTable_ref, $chr, @{$sorted_list}[$j], @{$sorted_list}[$k]) < $threshold) {
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

sub printUsage{
	print " Usage: \n";
	print " ./cgLD_Analysis_haploInfo_permuteRsq_v2.pl [rsq_out_put_name] [cg list]< haploInfo.txt\n";
	exit 0;
}

main();
