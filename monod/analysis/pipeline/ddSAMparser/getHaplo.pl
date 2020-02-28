#!/usr/bin/perl -w
# this script reads in a SAM format file with fillmd field, and output the pileup per read


use strict;
use warnings;

my %cgTable;
my %haploTable;
my $coverage_1 = 0;
my $coverage_2 = 0;
my $err1 = 0.01;
my $err2 = 0.001;

my $min_coverage = 10;
my $qual_base = 33;
my $min_length = 3;
my $rejection_cutoff = 0.1;
my $logp_for_1 = 0;
my $logp_for_2 = 0;

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
	while(my $line = <STDIN>){
		chomp($line);
		next if($line !~ m/MD:Z/);
		my %alignment = getAlignmentInfo($line);
		getReference(\%alignment);
		getCGPos(\%alignment);
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
		#print join(":", @ordered_pos), "\n";
		next if(scalar(@ordered_pos) < $min_length);
		my $posterior_1 = 1;
		my $posterior_2 = 1;
		my $total_1 = 0;
		my $total_2 = 0;
		my $list_pos = $alignment{"chr"};
		my $list_calls;
		for(my $i = 0; $i < scalar(@ordered_pos); $i++){
			my $pos = $ordered_pos[$i];
			$list_pos = defined($list_pos) ? $list_pos . ":" . $pos : $pos;
			$list_calls = defined($list_calls) ? $list_calls.$indices{$pos} : $indices{$pos};
		}
		$haploTable{$list_pos}->{$list_calls}++;
		undef(@ordered_pos);
		undef(%indices);	
	}
	undef(%cgTable);
	foreach my $value(keys %haploTable){
		my @can = keys %{$haploTable{$value}};
		my @posInfo = split ":", $value;
		my ($chr, $start_pos, $end_pos) = ($posInfo[0], $posInfo[1], $posInfo[scalar(@posInfo)-1]);
		foreach my $candidate (@can){
			print "$chr:$start_pos-$end_pos", "\t", $candidate, "\t", $haploTable{$value}->{$candidate}, "\t", join(":", @posInfo[1..scalar(@posInfo)-1]), "\n";
		}
	}
}

sub readCGTable{
        open(IN, "$ARGV[0]") || die("Error opening $ARGV[0]\n");
        while(my $line = <IN>){
                chomp($line);
                my @tmp = split /\t/, $line;
                $tmp[0] =~ s/:W//;
                $cgTable{$tmp[0].":".$tmp[1]} = 1;
        }
        close(IN);
}

sub getReference{
	my $alignment = shift;
	my $cigar = $alignment->{"cigar"};
	my $query = $alignment->{"seq"};
	my $qual = $alignment->{"qual"};
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
	#print "Ref_match_seq: $ref_seq\n";
	#print "Ref_match_qual: $ref_qual\n";

	$alignment->{"ref_match_seq"} = $ref_seq;
	$alignment->{"ref_match_qual"} = $ref_qual;
	$alignment->{"ref_nm_i"} = $nm_i;
	$alignment->{"ref_nm_d"} = $nm_d;
	$alignment->{"tlen"} = length($ref_seq) - $nm_i if($ref_seq);
}

sub getCGPos{
	my $alignment = shift;
	my @mismatches;
	my $mdz = $alignment->{"mdz"};
	my $ref_seq = uc($alignment->{"ref_match_seq"});
	my $ref_qual = $alignment->{"ref_match_qual"};
	return @mismatches if(!$alignment->{"ref_match_seq"});
	return @mismatches if(!$alignment->{"ref_match_qual"});
	#print $ref_seq, "\n", $ref_qual, "\n";
	my $chr = $alignment->{"chr"};
	my $position = $alignment->{"pos"};
	my @tmp = split /(\d+)/, $mdz;
	my @bases = split "", $ref_seq;
	my @quals = split "", $ref_qual;
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
				if($alignment->{"orientation"} eq "W"){
					push(@{$alignment->{"cgPos"}}, "$actual_pos:$b") if($cgTable{"$chr:$actual_pos"});
				}elsif($alignment->{"orientation"} eq "C"){
					$actual_pos--;
					die("Unrecognizable base $b\n") if(!$rcTable{$b});
					my $rc_b = $rcTable{$b};
					push(@{$alignment->{"cgPos"}}, "$actual_pos:$rc_b") if($cgTable{"$chr:$actual_pos"});
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
			if($alignment->{"orientation"} eq "W"){
				push(@{$alignment->{"cgPos"}}, "$actual_pos:$b") if($cgTable{"$chr:$actual_pos"});
			}elsif($alignment->{"orientation"} eq "C"){
				$actual_pos--;	
				die("Unrecognizable base $b\n") if(!$rcTable{$b});
				my $rc_b = $rcTable{$b};
				push(@{$alignment->{"cgPos"}}, "$actual_pos:$rc_b") if($cgTable{"$chr:$actual_pos"});
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
	$alignment{"ref_match_seq"} = "NA";
	$alignment{"ref_match_qual"} = "NA";
	$alignment{"ref_nm_i"} = 0;
	$alignment{"ref_nm_d"} = 0;
	$alignment{"tlen"} = 0;
	$alignment{"cgPos"} = ();
	return %alignment;
}

sub printUsage{
	print " Usage: \n";
	print " ./getHaplo.pl [cpg position list] < fillmd.sam\n";
	exit 0;
}

main();
