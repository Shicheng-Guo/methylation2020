#!/usr/bin/perl -w

use strict;
use List::Util qw(shuffle);
use Statistics::Descriptive;
use Statistics::R;

die("Not enough ARG") if(!$ARGV[3]);

my $scripts_dir = "/home/ddiep/scripts/ddSAMparser";
my $script;
$script = "$scripts_dir/mod.findClosestAndClusterHaploMaxLL.pl" if($ARGV[3] eq 'closest');
$script = "$scripts_dir/clusterHaploMaxLL.pl" if($ARGV[3] eq 'probeID');

my $ref_file_B = $ARGV[0];#cancer file
my $test_file = $ARGV[1]; #to create cross validation splits
my $num_folds = $ARGV[2];
my $use_method = $ARGV[4];
my $total_haps_in_test = 0;
my $minHaplo = 50;
my %badRegions;
my %mixTable;

sub main{
	my $sstime = time;
	my %splitFiles = generate_leave_one_out($test_file);
	my $entime = time;
	print "Elapsed time for generating leave one out: ", $entime-$sstime, "\n";
	$sstime = time;
	my %resultsTable = calculateAssignments(0, \%splitFiles);
	$entime = time;
	print "Elapsed time for calculate Assignment: ", $entime-$sstime, "\n";
	$sstime = time;
	my $stat_assigned = Statistics::Descriptive::Full->new();
	my $stat_correct = Statistics::Descriptive::Full->new();

	my @data_incorrect;
	my @data_assigned;
	# assign bad regions
	foreach my $region(keys %resultsTable){
		$resultsTable{$region}->{"HapRefB"} = 0 if(!$resultsTable{$region}->{"HapRefB"});
		$resultsTable{$region}->{"Unassigned"} = 0 if(!$resultsTable{$region}->{"Unassigned"});
		$resultsTable{$region}->{"HapUsed"} = 0 if(!$resultsTable{$region}->{"HapUsed"});
		push(@data_incorrect, $resultsTable{$region}->{"HapRefB"}/$resultsTable{$region}->{"HapUsed"});
		push(@data_assigned, 1-$resultsTable{$region}->{"Unassigned"}/$resultsTable{$region}->{"HapUsed"});
		if($resultsTable{$region}->{"HapRefB"}/$resultsTable{$region}->{"HapUsed"} > 0.25 ||
			1 - $resultsTable{$region}->{"Unassigned"}/$resultsTable{$region}->{"HapUsed"} < 0.5){	
			$badRegions{$region}=1;
		}
	}
	plotHistogram(\@data_incorrect, "Minor component rate");
	plotHistogram(\@data_assigned, "Assigned rate");
	$stat_assigned->add_data(@data_assigned);
	$stat_correct->add_data(@data_incorrect);
        print "Per region correct assigned mean: ", sprintf("%6.5f", $stat_correct->mean()), " stdev: ", sprintf("%6.5f", $stat_correct->standard_deviation()), "\n";
        print "Per region assigned mean: ", sprintf("%6.5f", $stat_assigned->mean()), " stdev: ", sprintf("%6.5f", $stat_assigned->standard_deviation()), "\n";
	$entime = time;
	print "Elapsed time for genering per region graphs/info: ", $entime - $sstime, "\n";

	# regenerate leave-1-out:
	%splitFiles = generate_leave_one_out($test_file);
	
	# Testing mixture
	$sstime = time;
	print "Testing mixtures\n";
	print "HapsMixed Unassigned_ave Unassigned_sd AssignedA_ave AssignedA_sd AssignedB_ave AssignedB_sd ScoreA_ave ScoreA_sd ScoreB_ave ScoreB_sd ScoreDiff_ave ScoreDiff_sd NumRegions\n";
	%resultsTable = calculateAssignments(0.000, \%splitFiles);
	$entime = time;
	print "Time for 0% mixture: ", $entime - $sstime, "\n";
	$sstime = time;
	%resultsTable = calculateAssignments(0.005, \%splitFiles);
	$entime = time;
	print "Time for 0.5% mixture: ", $entime - $sstime, "\n";
	$sstime = time;
	%resultsTable = calculateAssignments(0.010, \%splitFiles);
	$entime = time;
	print "Time for 1% mixture: ", $entime - $sstime, "\n";
	$sstime = time;
	%resultsTable = calculateAssignments(0.020, \%splitFiles);
	$entime = time;
	print "Time for 2% mixture: ", $entime - $sstime, "\n";
	$sstime = time;
	%resultsTable = calculateAssignments(0.040, \%splitFiles);
	$entime = time;
	print "Time for 4% mixture: ", $entime - $sstime, "\n";
	$sstime = time;
	%resultsTable = calculateAssignments(0.080, \%splitFiles);
	$entime = time;
	print "Time for 8% mixture: ", $entime - $sstime, "\n";
	$sstime = time;
	%resultsTable = calculateAssignments(0.160, \%splitFiles);
	$entime = time;
	print "Time for 16% mixture: ", $entime - $sstime, "\n";

	open(OUT_DATA, ">TestMixturesOut.txt") || die("Error writing to TestMixturesOut.txt");
	print OUT_DATA "Mix\tFrac.Minor\tScoresMajor\tScoresMinor\n";
	foreach my $mix (keys %mixTable){
		my @assigned_B = @{$mixTable{$mix}->{"assigned_B"}};
		my @scores_A = @{$mixTable{$mix}->{"scores_A"}};
		my @scores_B = @{$mixTable{$mix}->{"scores_B"}};
		for(my $i = 0; $i < $num_folds; $i++){
			print OUT_DATA "$mix\t", $assigned_B[$i], "\t", $scores_A[$i], "\t", $scores_B[$i], "\n";
		}
	}
	close(OUT_DATA);
	
	undef(%splitFiles);
	undef(%resultsTable);
	undef(%badRegions);
}

sub calculateAssignments{
	my $mixture = shift;
	my $splitFiles = shift;
	my %resultsTable_tmp;
	my ($unassigned_sum, $unassigned_sum_sqr) = (0,0);
	my ($assigned_A_sum, $assigned_A_sum_sqr) = (0,0);
	my ($assigned_B_sum, $assigned_B_sum_sqr) = (0,0);
	my ($scores_A_sum, $scores_A_sum_sqr) = (0,0);
	my ($scores_B_sum, $scores_B_sum_sqr) = (0,0);
	my ($scores_diff_sum, $scores_diff_sum_sqr) = (0,0);

	undef($mixTable{$mixture}->{"unassigned"});
	undef($mixTable{$mixture}->{"assigned_A"});
	undef($mixTable{$mixture}->{"assigned_B"});
	undef($mixTable{$mixture}->{"scores_A"});
	undef($mixTable{$mixture}->{"scores_B"});
	undef($mixTable{$mixture}->{"scores_diff"});

	my $RefFile_B = $ref_file_B;
	my $num_sample_in_mix = int($mixture*$total_haps_in_test/$num_folds);
	if($num_sample_in_mix > 1){
		generate_split_files($ref_file_B, $num_sample_in_mix, "ref_B.$mixture",1);
		#$RefFile_B = "ref_B.$mixture.haploInfo.txt";
	}

	for(my $i = 0; $i < $num_folds; $i++){
		my $n = sprintf("%03d", $i);
		my $RefFile_A = "ref_A.haploInfo.txt";
		my $TestFile = "ref_A.piece.haploInfo.txt";
		open(OUT, ">$RefFile_A") || die("Error writing to $RefFile_A\n");
		print OUT join("\n", @{${$splitFiles}{$n}->{"keep"}});
		close(OUT);
		open(OUT, ">$TestFile") || die("Error writing to $TestFile\n");
		print OUT join("\n", @{${$splitFiles}{$n}->{"out"}});
		close(OUT);
		my $MixTestFile = $TestFile;
		if($num_sample_in_mix > 1){
			# generate mixtures
			$MixTestFile = "ref_A.mixedpiece.haploInfo.txt";
			generate_split_files($TestFile, $num_sample_in_mix, "test.$mixture", 0);
			system("cat test.$mixture.haploInfo.txt ref_B.$mixture.sample.haploInfo.txt > $MixTestFile");
			unlink("test.$mixture.haploInfo.txt");
		}elsif($mixture > 0){
			print "Will not test mix with 1 or fewer\n";
		}

		my $cmd = "$script $RefFile_A $RefFile_B $MixTestFile $use_method";
		my $results = `$cmd`;
		chomp($results);
		my @fields = split /\ /, $results;
		my $countA = $fields[3];
		my $countB = $fields[4];
		my $unassigned = $fields[5];
		my $scoreA = $fields[8];
		my $scoreB = $fields[11]; 
		$countA =~ s/RefA=//g;
		$countB =~ s/RefB=//g;
		$unassigned =~ s/Unassigned=//g;

		die("Error assigning reads.\n") if($countA + $countB == 0);
		my $rateA = $countA/($countA+$countB);
		my $rateB = $countB/($countA+$countB);
	
		$unassigned_sum += $unassigned;
		$assigned_A_sum += $rateA;
		$assigned_B_sum += $rateB;
		$scores_A_sum += $scoreA;
		$scores_B_sum += $scoreB;
		$scores_diff_sum += ($scoreB-$scoreA);
	
		push(@{$mixTable{$mixture}->{"unassigned"}}, $unassigned);
		push(@{$mixTable{$mixture}->{"assigned_A"}}, $rateA);
		push(@{$mixTable{$mixture}->{"assigned_B"}}, $rateB);
		push(@{$mixTable{$mixture}->{"scores_A"}}, $scoreA);
		push(@{$mixTable{$mixture}->{"scores_B"}}, $scoreB);
		push(@{$mixTable{$mixture}->{"scores_diff"}}, ($scoreB-$scoreA));
		
		$unassigned_sum_sqr += $unassigned*$unassigned;
		$assigned_A_sum_sqr += $rateA*$rateA;
		$assigned_B_sum_sqr += $rateB*$rateB;
		$scores_A_sum_sqr += $scoreA*$scoreA;
		$scores_B_sum_sqr += $scoreB*$scoreB;
		$scores_diff_sum_sqr += ($scoreB-$scoreA)*($scoreB-$scoreA);

		open(IN_DATA, "MinorHaps.txt") || die("Error reading MinorHaps.txt\n");
		while(my $hapLine = <IN_DATA>){
			my ($probeID, $countB, $unassigned, $total) = split /\t/, $hapLine;
			$resultsTable_tmp{$probeID}->{"HapUsed"}+=$total;
			$resultsTable_tmp{$probeID}->{"HapRefB"}+=$countB;
			$resultsTable_tmp{$probeID}->{"Unassigned"}+=$unassigned;
		}
		close(IN_DATA);
		#print "My results: $results\n";		
	}
	unlink("ref_B.$mixture.haploInfo.txt");
	unlink("ref_B.$mixture.sample.haploInfo.txt");
	unlink("ref_A.mixedpiece.haploInfo.txt");
	unlink("ref_A.piece.haploInfo.txt");
	unlink("ref_A.haploInfo.txt");
	
	print "$num_sample_in_mix";
	print " ", sprintf("%.2f", $unassigned_sum/$num_folds), " ", sprintf("%.2f", naive_stdev($unassigned_sum, $unassigned_sum_sqr, $num_folds)),  
	      " ", sprintf("%.2f", $assigned_A_sum/$num_folds), " ", sprintf("%.2f", naive_stdev($assigned_A_sum, $assigned_A_sum_sqr, $num_folds)),
	      " ", sprintf("%.2f", $assigned_B_sum/$num_folds), " ", sprintf("%.2f", naive_stdev($assigned_B_sum, $assigned_B_sum_sqr, $num_folds)),
	      " ", sprintf("%.2f", $scores_A_sum/$num_folds), " ", sprintf("%.2f", naive_stdev($scores_A_sum, $scores_A_sum_sqr, $num_folds)),
	      " ", sprintf("%.2f", $scores_B_sum/$num_folds), " ", sprintf("%.2f", naive_stdev($scores_B_sum, $scores_B_sum_sqr, $num_folds)),
	      " ", sprintf("%.2f", $scores_diff_sum/$num_folds), " ", sprintf("%.2f", naive_stdev($scores_diff_sum, $scores_diff_sum_sqr, $num_folds)),
	      " ", scalar(keys %resultsTable_tmp), "\n";
	return %resultsTable_tmp;
}

sub naive_stdev{
	my $sum = shift;
	my $sum_sqr = shift;
	my $n = shift;
	return sqrt($sum_sqr - ($sum*$sum)/$n)/sqrt($n-1);
}

sub plotHistogram{
	my $data = shift;
	my $title = shift;
	my $output_file = $title.".png";
	$output_file =~ s/\ /_/g;
	my $R = Statistics::R->new();
	$R->startR;
	$R->set('x', $data);
	$R->set('chart_title', $title);
	$R->run(qq`png("$output_file")`);
	$R->run(q`hist(x, 100, main=chart_title, xlab="Rate")`);
	$R->run(q`dev.off()`);
	$R->stop();
}

#make mixture files
sub generate_split_files{
	my $file = shift;
	my $num_haplo_in_sample = shift;
	my $prefix = shift;
	my $take_sample = shift;
	my @numbers;
	my @lines_in_file;
	my ($line, $n, $i);
	#print "Splitting $file\n";
	open(OUT1, ">$prefix.haploInfo.txt") || die("Error writing to file\n");
	if($take_sample){
		open(OUT2, ">$prefix.sample.haploInfo.txt") || die("Error writing to file\n");
	}
	open(IN, "$file") || die("Error reading $file\n");
	$i = 0;
	while($line = <IN>){
		chomp($line);
		my ($probeID, $hapString, $count, $cpgPos) = split /\t/, $line;
		next if($badRegions{$probeID});
		$n = 0;
		push(@lines_in_file, $line);
		do{
			push(@numbers, $i);
			$i++;
			$n++;
		}while($n < $count);
	}
	close(IN);
	#print "\tTotal haps in $file: $i \n";
	my @shuffledNumbers = shuffle(@numbers);
	my @keep = sort {$a<=>$b} @shuffledNumbers[0..$num_haplo_in_sample-1];
	$i = 0;
	my $j = shift(@keep);
	$line = shift(@lines_in_file);
	do{
		my $keep_i = 0;
		$n = 0;
		my ($probeID, $hapString, $count, $cpgPos) = split /\t/, $line;
		do{
			if($j and $i == $j){
				$keep_i++; 
				$j = shift(@keep);
			}
			$i++;
			$n++;
		}while($n < $count);
		print OUT1 "$probeID\t$hapString\t", $count-$keep_i, "\t$cpgPos\n" if($count - $keep_i > 0);
		print OUT2 "$probeID\t$hapString\t$keep_i\t$cpgPos\n" if($keep_i > 0  and $take_sample);
		$line = shift(@lines_in_file);
	}while($line);
	close(OUT1);
	close(OUT2) if($take_sample);
}

#make leave-1-out file
sub generate_leave_one_out{
	my $file = shift;
	my %splitFiles;
	my @numbers;
	$total_haps_in_test = 0;
	for(my $i = 0; $i < $num_folds; $i++){
		my $n = sprintf("%03d", $i);
		push(@numbers, $n);
	}
	open(IN, "$file") || die("Error reading $file\n");
	while(my $line = <IN>){
		chomp($line);
		my ($probeID, $hapString, $count, $cpgPos) = split /\t/, $line;
		next if($badRegions{$probeID});
		$total_haps_in_test+=$count;
		my %keep;
		my %leftout;
		my $i = $count;
		do{
			my @shuffledNumbers = shuffle(@numbers);
			my $n = pop(@shuffledNumbers);
			$leftout{$n}++;
			foreach my $not_n (@shuffledNumbers){
				$keep{$not_n}++;
			}
			$i--;
		}while($i > 0);
		foreach my $n (keys %keep){
			push(@{$splitFiles{$n}->{"keep"}}, "$probeID\t$hapString\t".$keep{$n}."\t$cpgPos");
		}
		foreach my $n (keys %leftout){
			push(@{$splitFiles{$n}->{"out"}}, "$probeID\t$hapString\t".$leftout{$n}."\t$cpgPos");
		}
	}
	return %splitFiles;
}
main();
