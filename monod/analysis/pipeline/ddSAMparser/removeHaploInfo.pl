#!/usr/bin/perl -w

use strict;
use List::Util qw(shuffle);

my $remove_num = $ARGV[0];
my $total_in_file = $ARGV[1];
die("Removing too many!!") if($total_in_file < $remove_num);

my @numbers;
push(@numbers, $_) for (1..$total_in_file);
my @shuffledNumbers = shuffle(@numbers);
my @skip = @shuffledNumbers[0..$remove_num-1];
#print scalar(@skip), "\n";
#print join("\n", @skip), "\n";
my @sorted_skip = sort{$a<=>$b} @skip;
my $cur_skip = 0;
#chr22:45103758-45103899 CCCCCC  1       45103765,45103814,45103836,45103846,45103849,45103857
my $tracker = 0;
while(my $line = <STDIN>){
	chomp($line);
	my ($probeID, $hapString, $count, $cpgPos) = split /\t/, $line;
	my $i = $count;
	my $keep = 0;
	do{
		$tracker++;
		if($tracker == $sorted_skip[$cur_skip]){
			#don't increment keep;
			$cur_skip++ if($cur_skip < $remove_num-1);
		}else{
			$keep++;
		}
		$i--;
	}while($i > 0);
	if($keep > 0){
		print "$probeID\t$hapString\t$keep\t$cpgPos\n";
	}
}

