#!/usr/bin/perl -w
# this script reads in a SAM format file with fillmd field, and output the pileup per read


use strict;
use warnings;

my %alignment;
my $qual_base = 33;
my $deletion_string = "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD";
my $padding_string = "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP";

sub main{
	printUsage() if(!$ARGV[0]);
	open(OUT, ">$ARGV[0]") || croak("error writing to file $ARGV[0]\n");
	while(my $line = <STDIN>){
		chomp($line);
		next if($line !~ m/MD:Z/);
		%alignment = getAlignmentInfo($line);
		getReference();
		my @mismatches = getMismatchPos();
		print OUT "numer of indel: i=", $alignment{"ref_nm_i"}, ", d=", $alignment{"ref_nm_d"}, "\n";
		print OUT "number of mm: ", scalar(@mismatches), "\n";
		foreach my $mm (@mismatches){
			print OUT "\t", $mm, "\n";
		}
		undef(%alignment);
		
	}
	close(OUT);
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
		die("Could not understand $op: '$cigar'") if($op !~ m/[MX=HDISP]/);
		$op =~ s/[X=]/M/g;
		my ($clip_s, $clip_q);
		if($op eq "M" || $op eq "I" || $op eq "S"){
			$clip_s = substr($query, $j, $runlen);
			$clip_q = substr($qual, $j, $runlen);
			$clip_s =~ s/[ATGCatgc]/I/g if($op eq "I");
			$nm_i += $runlen if($op eq "I");
			$j += $runlen;
		}elsif($op eq "H"){
			#do nothing;
		}elsif($op eq "P"){
			$clip_s = substr($padding_string, 0, $runlen);
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
}

sub getMismatchPos{
	my $mdz = $alignment{"mdz"};
	my $ref_seq = $alignment{"ref_match_seq"};
	my $ref_qual = $alignment{"ref_match_qual"};
	my $position = $alignment{"pos"};
	my @tmp = split /(\d+)/, $mdz;
	my @bases = split "", $ref_seq;
	my @quals = split "", $ref_qual;
	my @mismatches;
	my $j = 0;
	my $mm_pos = 0;
	for(my $i = 1; $i< scalar(@tmp); $i++){
		my $op = $tmp[$i];
		if($op =~ s/^\^//){
			$j+=length($op);
			$mm_pos+=length($op);
		}elsif($op =~ /^[+-]?\d+$/){
			for(my $k = 0; $k < $op; $k++){
				$mm_pos++ if($bases[$j] ne "I");
			}
			$j+=$op;
		}elsif($op =~ /[ATGCatgc]/){
			my $b = $bases[$j];
			my $q = $quals[$j];
			my $actual_pos = $position + $mm_pos;
			my $score = ord($q) - $qual_base;
			push(@mismatches, "$b, $op, $score, $actual_pos");
			$mm_pos++;
			$j++;
		}
	}
	return @mismatches;
}

sub getAlignmentInfo{
	my $line = shift;
	my @tmp = split /\t/, $line;
	my %alignment;
	$alignment{"name"} = $tmp[0];
	$alignment{"flag"} = $tmp[1];
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
	print " ./ddSAMparse.pl [out_put_name] < fillmd.sam\n";
	exit 0;
}

main();
