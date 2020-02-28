#!/usr/bin/perl -w
# this script reads in a SAM format file with fillmd field, and output the pileup per read


use strict;
use warnings;

my %cgTable;
my %methylTable;
my $coverage_1 = 0;
my $coverage_2 = 0;
my $err1 = 0.01;
my $err2 = 0.001;

my $min_coverage = 10;
my $qual_base = 33;
my $min_length = 5;
my $rejection_cutoff = 0.2;
my $logp_for_1 = 0;
my $logp_for_2 = 0;
my $file1_count = 0;
my $file2_count = 0;

my %alignment;
my $deletion_string = "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD";

my %rcTable;
$rcTable{'A'}='T';
$rcTable{'T'}='A';
$rcTable{'G'}='C';
$rcTable{'C'}='G';
$rcTable{'N'}='N';
$rcTable{'D'}='D';

my @palette=("0,240,0", "30,210,0","60,180,0","90,150,0","120,120,0","150,90,0","180,60,0","210,0,0");

sub main{
	printUsage() if(!$ARGV[0]);
	readBED_CGTable();
	#print "ReadID\tIndices\tAlleleCalls\tPosterior_1\tPosterior_2\n";
	while(my $line = <STDIN>){
		chomp($line);
		next if($line !~ m/MD:Z/);
		%alignment = getAlignmentInfo($line);
		getReference();
		getCGPos();
		next if(!$alignment{"cgPos"});
		my @cgPos = @{$alignment{"cgPos"}};
		my %indices;
		foreach my $pos_call (@cgPos){
			my ($pos, $call) = split ":", $pos_call;
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
		next if(scalar(@ordered_pos) < $min_length);
		my $posterior_1 = 1;
		my $posterior_2 = 1;
		my $total_1 = 0;
		my $total_2 = 0;
		my $list_pos = $alignment{"chr"};
		my $list_calls;
		for(my $i = 0; $i < scalar(@ordered_pos); $i++){
			my $pos = $ordered_pos[$i];
			my $index = $alignment{"chr"}.":".$pos;
			next if(!$cgTable{$index});
			my ($f1U, $f1M, $f2U, $f2M) = ($methylTable{$index}->{"1U"}, 
							$methylTable{$index}->{"1M"},
							$methylTable{$index}->{"2U"},
							$methylTable{$index}->{"2M"});
			my $total = $f1U + $f1M + $f2U + $f2M;
			$total_1 += $f1U + $f1M;
			$total_2 += $f2U + $f2M;
			my ($p1U, $p1M, $p2U, $p2M, $pU, $pM) = ($f1U/($f1U+$f1M), 
									$f1M/($f1U+$f1M),
									$f2U/($f2U+$f2M),
									$f2M/($f2U+$f2M),
									($f1U+$f2U)/$total,
									($f1M+$f2M)/$total);
			if($indices{$pos} eq "C"){
				$posterior_1 *= ((1-$err1)*$p1M + $err1*$p1U);
				$posterior_2 *= ((1-$err1)*$p2M + $err1*$p2U);
			}else{
				$posterior_1 *= ((1-$err2)*$p1U + $err2*$p1M);
				$posterior_2 *= ((1-$err2)*$p2U + $err2*$p2M);
			}
			$list_pos = $list_pos . ":" . $pos;
			$list_calls = defined($list_calls) ? $list_calls.$indices{$pos} : $indices{$pos};
		}
		if(defined($list_calls) and length($list_calls) >= $min_length){
			$posterior_1 *= $total_1/($total_1+$total_2);
			$posterior_2 *= $total_2/($total_1+$total_2);
			#next if($posterior_1 < $rejection_cutoff and $posterior_2 < $rejection_cutoff);
			#$posterior_1 *= 0.5;
			#$posterior_2 *= 0.5;
			my $evidence = $posterior_1+$posterior_2;
			if($evidence and $evidence > $rejection_cutoff){
				$posterior_1 /= $evidence;
				$posterior_2 /= $evidence;
				$file1_count++ if($posterior_1 > $posterior_2);
				$file2_count++ if($posterior_2 > $posterior_1);
				#print $alignment{"name"}, "\t", $list_pos, "\t", $list_calls, "\t", $posterior_1,"\t", $posterior_2, "\t", $evidence, "\n";
				$logp_for_1 += log($posterior_1);
				$logp_for_2 += log($posterior_2);
			}
		}
		undef(@ordered_pos);
		undef(%indices);
		undef(%alignment);	
	}
	print "File 1 posterior=$logp_for_1\n";
	print "File 1 reads=$file1_count\n";
	print "File 2 posterior=$logp_for_2\n";
	print "File 2 reads=$file2_count\n";
	undef(%cgTable);
	undef(%methylTable);
}

sub readBED_CGTable{
	# cgTable{chr:pos}
	# methylTable{chr:pos}
	#chr10	100017297	100017298	'194/200'	200	+	100017297	100017298	210,0,0

	open(IN, "$ARGV[0]") || die("Error opening $ARGV[0]\n");
	while(my $line = <IN>){
		next if($line =~ /track/);
		chomp($line);
		my @f = split /\t/, $line;
		my $index = $f[0].":".$f[2];
		$f[3] =~ s/'//g;
		my ($c, $depth) = split /[\/]/, $f[3];
		next if($depth < $min_coverage);
		$methylTable{$index}->{"1M"} = $c;
		$methylTable{$index}->{"1U"} = $depth - $c;
		$coverage_1+=$depth;
	}
	close(IN);

	open(IN, "$ARGV[1]") || die("Error opening $ARGV[1]\n");
	while(my $line = <IN>){
		next if($line =~ /track/);
		chomp($line);
		my @f = split /\t/, $line;
		my $index = $f[0].":".$f[2];
		next if(!$methylTable{$index});
		$f[3] =~ s/'//g;
		my ($c, $depth) = split /[\/]/, $f[3];
		next if($depth < $min_coverage);
		$cgTable{$index}=1;
		$methylTable{$index}->{"2M"} = $c;
		$methylTable{$index}->{"2U"} = $depth - $c;
		$coverage_2+=$depth;
	}
	close(IN);
}

sub getReference{
	my $cigar = $alignment{"cigar"};
	my $query = $alignment{"seq"};
	my $qual = $alignment{"qual"};
	my $ref_seq;
	my $ref_qual;
	# Adapted from Ben Langmead (BtlBio::Alignment:Util.pm)
	# CIGAR fields are always in pairs
	my $i = 0;
	my $j = 0;
	my $nm_i = 0;
	my $nm_d = 0;
	while($i < length($cigar)){
		substr($cigar, $i) =~ /^([0-9]+)/;
		defined($1) || die("Could not parse number at pos $i: '$cigar'");
		my $runlen = $1;
		$i += length($1);
		$i < length($cigar) || die("Bad cigar string : '$cigar'");
		my $op = substr($cigar, $i, 1);
		defined($op) || die("count not parse operation at pos $i: '$cigar'");
		die("Could not understand $op: '$cigar'") if($op !~ m/[MX=DIS]/);
		$op =~ s/[X=]/M/g;
		my ($clip_s, $clip_q);
		if($op eq "M" || $op eq "I" || $op eq "S"){
			$clip_s = substr($query, $j, $runlen);
			$clip_q = substr($qual, $j, $runlen);
			$clip_s =~ s/[ATGCatgc]/I/g if($op eq "I");
			$nm_i += $runlen if($op eq "I");
			$j += $runlen;
		}else{
			#deletion from reference
			$nm_d += $runlen;
			length($deletion_string) > $runlen || die("deletion is too long at $runlen: '$cigar'");
			$clip_s = substr($deletion_string, 0, $runlen);
			$clip_q = substr($deletion_string, 0, $runlen);
		}
		$i++;
		$ref_seq = $ref_seq . $clip_s if($op =~ m/[MD]/);
		$ref_qual = $ref_qual . $clip_q if($op =~ m/[MD]/);
	}
	$alignment{"ref_match_seq"} = $ref_seq;
	$alignment{"ref_match_qual"} = $ref_qual;
	$alignment{"ref_nm_i"} = $nm_i;
	$alignment{"ref_nm_d"} = $nm_d;
	$alignment{"tlen"} = length($ref_seq) - $nm_i;
}

sub getCGPos{
	my $mdz = $alignment{"mdz"};
	my $ref_seq = uc($alignment{"ref_match_seq"});
	my $ref_qual = $alignment{"ref_match_qual"};
	my $chr = $alignment{"chr"};
	my $position = $alignment{"pos"};
	my @tmp = split /(\d+)/, $mdz;
	my @bases = split "", $ref_seq;
	my @quals = split "", $ref_qual;
	my @mismatches;
	my $j = 0;
	my $relative_pos = 0;
	for(my $i = 1; $i< scalar(@tmp); $i++){
		my $op = $tmp[$i];
		if($op =~ s/^\^//){
			$j+=length($op);
			$relative_pos+=length($op);
			#print $op, "\n";
		}elsif($op =~ /^[+-]?\d+$/){
			for(my $k = 0; $k < $op; $k++){
				# obtain calls at the current position in ref_seq
				my $b = $bases[$j];
				my $q = $quals[$j];
				my $actual_pos = $position + $relative_pos;
				$j++; # j position moves along ref_seq
				if($b eq "I"){ $k--; next;} # if the position is an insertion to the reference, we need to go back 1
				$relative_pos++; # relative_position on get incremented when the position is not an insertion to the reference
				next if(ord($q) - $qual_base < 5); # do not make a call if the quality is bad
				#print "SM: $relative_pos\t$op\t$b\t$q\n";
				if($alignment{"orientation"} eq "W"){
					push(@{$alignment{"cgPos"}}, "$actual_pos:$b") if($cgTable{"$chr:$actual_pos"});
				}elsif($alignment{"orientation"} eq "C"){
					$actual_pos--;
					die("Unrecognizable base $b\n") if(!$rcTable{$b});
					my $rc_b = $rcTable{$b};
					push(@{$alignment{"cgPos"}}, "$actual_pos:$rc_b") if($cgTable{"$chr:$actual_pos"});
				}
			}
		}elsif($op =~ /[ATGCatgc]/){
			# obtain calls at the current position in ref_seq
			my $b = $bases[$j];
			my $q = $quals[$j];
			my $actual_pos = $position + $relative_pos;
			# make the necessary incrementations in j and relative positions
			$relative_pos++;
			$j++;
			next if(ord($q) - $qual_base < 5); # do not make a call if the quality is bad
			#print "NM: $relative_pos\t$op\t$b\t$q\n";
			if($alignment{"orientation"} eq "W"){
				push(@{$alignment{"cgPos"}}, "$actual_pos:$b") if($cgTable{"$chr:$actual_pos"});
			}elsif($alignment{"orientation"} eq "C"){
				$actual_pos--;	
				die("Unrecognizable base $b\n") if(!$rcTable{$b});
				my $rc_b = $rcTable{$b};
				push(@{$alignment{"cgPos"}}, "$actual_pos:$rc_b") if($cgTable{"$chr:$actual_pos"});
			}
		}
	}
	undef @bases;
	undef @quals;
	undef @tmp;
	return @mismatches;
}

sub getAlignmentInfo{
	my $line = shift;
	my @tmp = split /\t/, $line;
	my %alignment;
	$alignment{"name"} = $tmp[0];
	$alignment{"orientation"} = $tmp[1] & 16 ? "C":"W";
	$alignment{"chr"} = $tmp[2];
	$alignment{"pos"} = $tmp[3];
	$alignment{"cigar"} = $tmp[5];
	$alignment{"seq"} = $tmp[9];
	$alignment{"qual"} = $tmp[10];
	$alignment{"ndz"} = $tmp[11];
	$alignment{"mdz"} = $tmp[12];
	return %alignment;
}

sub printUsage{
	print " Usage: \n";
	print " ./naiveBayes_cgBED.pl [BED1] [BED2] < fillmd.sam\n";
	exit 0;
}

main();
