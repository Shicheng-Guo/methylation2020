#!/usr/bin/perl -w
# this script reads in a SAM format file with fillmd field, and output the pileup per read


use strict;
use warnings;

my %cgTable;
my %allelesTable;
my @alleles;
my $minCoverage = 10;
my $qual_base = 33;

my $locus_num_cg;

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
	readCGTable();
	readCGAllelesTable();
	open(OUT, ">$ARGV[0]") || die("error writing to file $ARGV[0]\n");
	print OUT "ReadID\tIndices\tAlleleCalls\tProbabilities\tTotalProbabilities\n";
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
		next if(scalar(@ordered_pos) < $locus_num_cg);
		my $total_probs = 1;
		my $num_locus = 0;
		my ($list_indices, $list_alleles, $list_probs);
		for(my $i = 0; $i <= scalar(@ordered_pos)-$locus_num_cg; $i++){
			my $value_key = $alignment{"chr"};
			my $value_call;
			for(my $j = 0; $j < $locus_num_cg; $j++){
				$value_key = $value_key . ":" . $ordered_pos[$i+$j];
				$value_call = $value_call . $indices{$ordered_pos[$i+$j]};
			}
			$value_call =~ tr/T/U/;
			$value_call =~ tr/C/M/;
			if($allelesTable{$value_key}->{"total_num"}){
				my $total_evidence = $allelesTable{$value_key}->{"total_num"};
				my $num = $allelesTable{$value_key}->{$value_call};
				$num = 0 if(!$num);
				my $p = sprintf("%4.3f", $num/$total_evidence);
				$list_indices = defined($list_indices) ? $list_indices.",".$value_key : $value_key;
				$list_alleles = defined($list_alleles) ? $list_alleles.",".$value_call : $value_call;
				$list_probs = defined($list_probs) ? $list_probs.",".$p : $p;
				$total_probs *= $p;
				$num_locus++;
			}
		}
		if(defined($list_indices)){
			print OUT $alignment{"name"}, "\t", $list_indices, "\t", $list_alleles, "\t", $list_probs, "\t", $total_probs, "\n";
		}
		undef(@ordered_pos);
		undef(%indices);
		undef(%alignment);	
	}
	
	close(OUT);
	undef(%allelesTable);
	undef(%cgTable);
}

sub readCGAllelesTable{
	open(IN, "$ARGV[2]") || die("Error opening $ARGV[2]\n");
	while(my $line = <IN>){
		next if($line =~ /indice/);
		chomp($line);
		my @tmp = split /\t/, $line;
		my $index = $tmp[0];	
		my @counts = split ",", $tmp[1];
		my $total_depth = $tmp[2];
		$allelesTable{$index}->{"total_num"} = $total_depth;
		foreach my $value (@counts){
			next if(!$value || $value !~ /[UM]/);
			my ($spacer, $d, $a) = split /(\d+)/, $value;
			if(!$locus_num_cg){
				$locus_num_cg = length($a);
				generateAlleles();
			}
			$allelesTable{$index}->{$a} = $d;
		}
	}
	close(IN);
}

sub readCGTable{
	open(IN, "$ARGV[1]") || die("Error opening $ARGV[1]\n");
	while(my $line = <IN>){
		chomp($line);
		my @tmp = split /\t/, $line;
		$tmp[0] =~ s/:W//;
		$cgTable{$tmp[0].":".$tmp[1]} = 1;		
	}
	close(IN);
}

sub generateAlleles{
        @alleles = ("U", "M");
        my $i = 1;
        while($i < $locus_num_cg){
                my @new_alleles;
                foreach my $value (@alleles){
                        push(@new_alleles, $value ."U");
                        push(@new_alleles, $value ."M");
                }
                undef @alleles;
                @alleles = @new_alleles;
                undef @new_alleles;
                $i++;
        }
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
	print " ./compare_cgLocusAlleles.pl [cpg position list] [cgLocusAlleles]  < fillmd.sam\n";
	exit 0;
}

main();
