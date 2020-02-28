#!/usr/bin/perl -w
#This script takes sorted sam file from the STDIN, and remove redundant reads based on the mapping positions and UMIs.
# keep the longest

my %umi_table;
my $last_chr_pos=0;
my $n_count = 0;

my %umi_unique;
my %posTable;

my $sum_umi_unique=0;
my $m_value = 0;

while(my $line = <STDIN>){
	if($line =~ /^\@/){
		print $line;
		next;
	}
	chomp($line);
	my @fields = split(/\t/, $line);
	my $read_id = $fields[0];
	my $chr_pos = $fields[2] . "-". $fields[3];
	my $score = $fields[4];
	$read_id =~ s/-1/\_/g;
	my @parts = split(/\_+/, $read_id);
	my $umi_tag;
	foreach my $value(@parts){
		next if($value =~ m/[0-9]/);
		next if(length($value) < 6);
		$umi_tag = $value;
		$m_value = 4**length($value) if($m_value eq 0);
		$value =~ s/:F//;
		$value =~ s/:R//;
		$umi_unique{$value}++;
		$sum_umi_unique++;
	}
	die("Can't find $umi_tag\n") if(!$umi_tag);
        if($chr_pos ne $last_chr_pos && $last_chr_pos ne 0){
		my $k_count = scalar(keys %umi_table);
		#=== select non-clonal reads here
		# calculate my true n-value
		my $max_keep = 1;
		$n_value = int(-1 * $m_value * log(1 - $k_count/$m_value));
		$n_value = $n_count if($n_count < $n_value);
		my $num_keep = $max_keep*$k_count;
		my $incr = $max_keep;
		while($num_keep < $n_value){
			$incr++;
			foreach my $umi (keys %umi_table){
				my @cand = @{$umi_table{$umi}};
				$num_keep++ if(scalar(@cand) >= $incr);
				last if($num_keep >= $n_value);
			}
		}
		$max_keep = $incr;
		foreach my $umi(keys %umi_table){
			my @cand = @{$umi_table{$umi}};
			if(scalar(@cand) <= $max_keep){
			# print everything
				print join("\n", @cand), "\n";
			}else{
			# print only the best-scoring alignments
				my %score_hash;
				foreach my $can(@cand){
					my @tmp = split /\t/, $can;
					$score_hash{$can} = $tmp[4];
				}
				my $k = 0;
				foreach my $can (sort {$score_hash{$a} <=> $score_hash{$b}} keys %score_hash){
					print $can, "\n";
					$k++;
					last if($k == $max_keep);
				}
				undef %score_hash;
			}
		}
		#================================
                undef %umi_table;
                push(@{$umi_table{$umi_tag}}, $line);
                $last_chr_pos = $chr_pos;
		$n_count = 1;
		next;
        }
	if(!$umi_table{$umi_tag}){
		push(@{$umi_table{$umi_tag}}, $line);
		$last_chr_pos = $chr_pos;
		$n_count++;
		next;
	}
	my $last_line = $line;
	my @last_fields = split /\t/, $last_line;
	push(@{$umi_table{$umi_tag}}, $line);
	$last_chr_pos = $chr_pos;
	$n_count++;
}

my $k_count = scalar(keys %umi_table);
#=== select non-clonal reads here
my $max_keep = 1;
$n_value = int(-1 * $m_value * log(1 - $k_count/$m_value));
$n_value = $n_count if($n_count < $n_value);
my $num_keep = $max_keep*$k_count;
my $incr = $max_keep;
while($num_keep < $n_value){
	$incr++;
	foreach my $umi (keys %umi_table){
		my @cand = @{$umi_table{$umi}};
		$num_keep++ if(scalar(@cand) >= $incr);
		last if($num_keep >= $n_value);
	}
}
$max_keep = $incr;
foreach my $umi(keys %umi_table){
	my @cand = @{$umi_table{$umi}};
	if(scalar(@cand) <= $max_keep){
	# print everything
		print join("\n", @cand), "\n";
	}else{
	# print only the best-scoring alignments
		my %score_hash;
		foreach my $can(@cand){
			my @tmp = split /\t/, $can;
			$score_hash{$can} = $tmp[4];
		}
		my $k = 0;
		foreach my $can (sort {$score_hash{$a} <=> $score_hash{$b}} keys %score_hash){
			print $can, "\n";
			$k++;
			last if($k == $max_keep);
		}
		undef %score_hash;
	}
}
undef %umi_table;
	
#=== print out the umi counts here
open(OUT, ">max_k_n_umi.txt") || die("Error writing to max_k_n.txt\n");
my $total_umi = scalar(keys %umi_unique);
foreach my $umi (keys %umi_unique){
	print OUT $umi, "\t", $umi_unique{$umi}, "\n";
}
close(OUT);
#================================
