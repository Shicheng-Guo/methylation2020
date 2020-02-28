#!/usr/bin/perl -w
#Usage:

# fixMateInfo.pl < SAM_SORTED_BY_NAME > STITCHED_SAM
# first, clean up the names and sam flags
# run bam clipOverlap
use strict;

#M01186:73:000000000-A20K5:1:1113:15735:6349_2:N:0:1_AAGTGT:F    129       chr10   1814747 60      124M    *       0       0       GAGTTATTGTGTTTGGTTAGGTGAATTTTATTGATGAGATGTTGTTGGGTTATTGGGATATTTGTAGTAGAAATGTTTATAAGGTAATAAAATGTGTAGTTTATGTAGTTTATGTTTTAGATTT    E00DD2FGGDFGFFCCGGHHHDG1AFGHHGHHHHHH2AAGHHHHHHFHGGFAGGHEH>F0FGHBGGHHGHF2BFGHHHFGHGDGHHHHHF1FFGHHHHHHHGHHHHHHHGHFHGHFHHHGHGGH

#M01186:73:000000000-A20K5:1:1113:15735:6349:1:N:0:1_AAGTGT:R    65       chr10   1814789 60      113M    *       0       0       TGTTGGGTTATTGGGATTTTTGTAGTAGAAATGTTTATAAGGTAATAAAATGGTTATTTTATTTAGTTTATGTTTTAGATTTTTTTTTGGGTTTGTTTTTTTTTTTTTTTGTT       HHAHGGHHGGEGGGGGGCEAGGFE/FGEHHE/@CA12BDGB1/1?GH1211?2?FF12111<111F11111<<1111?F11111111>11>11FB...//0<</<<0C.C..C

sub main{
	my %prev_alignment;
	while(my $line = <STDIN>){
		chomp($line);
		my %alignment = getAlignmentInfo($line);
		if($prev_alignment{"name"}){
			if($prev_alignment{"name"} eq $alignment{"name"} and 
			   $prev_alignment{"flag"} eq $alignment{"flag"} and 
			   $prev_alignment{"chr"} eq $alignment{"chr"}){
				$prev_alignment{"mate_chr"} = $alignment{"chr"};
				$prev_alignment{"mate_pos"} = $alignment{"pos"};
				$alignment{"mate_chr"} = $prev_alignment{"chr"};
				$alignment{"mate_pos"} = $prev_alignment{"pos"};
				if($prev_alignment{"flag"} eq 0){
					$prev_alignment{"flag"} = $prev_alignment{"pos"} <= $alignment{"pos"} ? 65 : 129;
					$alignment{"flag"} = $alignment{"pos"} < $prev_alignment{"pos"} ? 65 : 129;
				}else{
					$prev_alignment{"flag"} = $prev_alignment{"pos"} <= $alignment{"pos"} ? 113 : 177;
					$alignment{"flag"} = $alignment{"pos"} < $prev_alignment{"pos"} ? 113 : 177;
				}
				print getLine(\%prev_alignment), "\n";
				print getLine(\%alignment), "\n";
				for (keys %prev_alignment){
					delete $prev_alignment{$_};
				}
			}else{
				print getLine(\%prev_alignment), "\n";
				for (keys %alignment){
					$prev_alignment{$_} = $alignment{$_};
				}
			}
		}else{
			%prev_alignment = %alignment;
		}
	}
	# process the last alignment:
	print getLine(\%prev_alignment), "\n" if($prev_alignment{"name"});	
}

sub getAlignmentInfo(){
        my $line = shift;
        my @tmp = split /\t/, $line;
        my %alignment;
	$tmp[0] =~ s/\/1$//g;
	$tmp[0] =~ s/\/2$//g;
	$tmp[0] =~ s/\/3$//g;
	if($tmp[0] =~ /_/){
		my @f = split "_", $tmp[0];
        	$alignment{"name"} = $f[0]."_".$f[1];
	}else{
		$alignment{"name"} = $tmp[0];
	}	
        $alignment{"flag"} = $tmp[1];
        $alignment{"chr"} = $tmp[2];
        $alignment{"pos"} = $tmp[3];
	$alignment{"score"} = $tmp[4];
        $alignment{"cigar"} = $tmp[5];
	$alignment{"mate_chr"} = $tmp[6];
	$alignment{"mate_pos"} = $tmp[7];
	$alignment{"fraglen"} = $tmp[8];
        $alignment{"seq"} = $tmp[9];
        $alignment{"qual"} = $tmp[10];
        return %alignment;
}

sub getLine(){
	my $alignment = shift;
	my $line = ${$alignment}{"name"}."\t".${$alignment}{"flag"}."\t".${$alignment}{"chr"}."\t".${$alignment}{"pos"}."\t".${$alignment}{"score"}."\t".${$alignment}{"cigar"}."\t".${$alignment}{"mate_chr"}."\t".${$alignment}{"mate_pos"}."\t".${$alignment}{"fraglen"}."\t".${$alignment}{"seq"}."\t".${$alignment}{"qual"};
	return $line;
}

main();
