#!/usr/bin/perl -w

# bisReadMapper.pl: a perl script to map bisulfite sequencing reads.
# This version performs soft trimming and uses SOAP/Bowtie2. Reads can either be trimmed & encoded already or not.
# Contact: Dinh Diep
# Version 1.3

use strict;
use warnings;
use Getopt::Std;

my @reads = ();
my ($encodedFqName1, $encodedFqName2);
my $align_mode = 'S';
my $cpu = 2;
my $fqName = "Sequence1";
my $trim3=0;
my $trim5=0;
my $qualtrim=0;
my $nomap = 0;
my $noencode = 0;
my $keep_bam = 1;
my $TMP_DIR;
my $qual_base = 64;
my $maxMismatches=5;
my $score_min= "L,-0.6,-0.6";
my $gem_allow="-m 0.04 -e 0.04 -s 0 --unique-mapping";
my $min_identity = 0.95;

my ($template_fwd, $template_rev, $template_idx);
my ($soap2_exe, $bowtie2_exe, $bwa_exe, $last_dir, $gem_dir) = (0,0,0,0,0);

my %rcTable;
$rcTable{'A'}='T';
$rcTable{'T'}='A';
$rcTable{'G'}='C';
$rcTable{'C'}='G';
$rcTable{'N'}='N';
$rcTable{'R'}='Y';
$rcTable{'Y'}='R';
$rcTable{'M'}='K';
$rcTable{'K'}='M';
$rcTable{'S'}='S';
$rcTable{'W'}='W';

my %chrSizes;
my %countReads;
my $missed_guesses = 0;
my $total_reads = 0;
my $total_bases = 0;
my $total_mbases = 0;

&main;
exit;

sub main(){
	my %opts = ();
	getopts('r:W:C:g:a:e:p:n:3:5:q:T:b:s:', \%opts);
	die( printUsage() ) if(@ARGV == 0 && !$opts{r});
	die("No reads provided.\n") if(!$opts{r});
	@reads = split ",", $opts{r};
	print "Reads: ", join(",", @reads), "\n";
	$fqName = $opts{n} if($opts{n});
	print "SampleName: ", $fqName, "\n";
	if(scalar(@reads) == 2){ $align_mode = 'P'; print "Aligning as paired-ends\n"};
	die("Too many reads files or commas in read file names not acceptable\n") if(scalar(@reads) > 2);
	if(!$opts{a}){ $nomap = 1; print "No aligner given\n"}
	$bowtie2_exe = $opts{a} if($opts{a} =~ m/bowtie/);
	$soap2_exe = $opts{a} if($opts{a} =~ m/soap/);
	$bwa_exe = $opts{a} if($opts{a} =~ m/bwa/);
	$last_dir = $opts{a} if($opts{a} =~ m/last/);
	$gem_dir = $opts{a} if($opts{a} =~ m/gem/i);
	if(!$nomap){ print "Aligning with ", $opts{a}, "\n"}
	$qual_base = $opts{b} if($opts{b});
	$score_min = $opts{s} if($opts{s});
	die("No reference index [Watson or Crick or both] given\n") if(!$opts{W} && !$opts{C});
	die("No reference genome fasta given\n") if(!$opts{g});
	$template_fwd = $opts{W};
	$template_rev = $opts{C} if($opts{C});
	$template_idx = $opts{g};
	$TMP_DIR = $opts{T};
	$cpu = $opts{p};
	if(!$TMP_DIR){
		my $cur_dir = `pwd`;
		chomp($cur_dir);
		$TMP_DIR = $cur_dir;
	}
	
	##--------------Encoding below-----------------##
	if($opts{e}){
		$noencode = 1;
		$encodedFqName1 = $reads[0];
		$encodedFqName2 = $reads[1] if($align_mode eq 'P');
	}else{
		$encodedFqName1 = $fqName . ".1.encoded";
		$encodedFqName2 = $fqName . ".2.encoded";
		$trim5 = $opts{5} if($opts{5});
		$trim3 = $opts{3} if($opts{3});
		$qualtrim = $opts{q} if($opts{q});
		encodeFastq($reads[0], $encodedFqName1, $trim5);
		encodeFastq($reads[1], $encodedFqName2, $trim5) if($align_mode eq 'P');
	}
	print "Reads: ", join(",", @reads), "\n";
	print "Finish encoding reads.\n";
	exit 0 if($nomap);

	##--------------Mapping below-----------------##
	die("Need to have faidx file: $template_idx\n") if(!getChromSizes());
	my $start_time = time;	
	my $map_file; # SAM formated map file.
	if($align_mode eq 'P'){
		$map_file = fastq2SOAPpe() if($soap2_exe);
		$map_file = fastq2BOWTIEpe() if($bowtie2_exe);
		$map_file = fastq2BWApe() if($bwa_exe);
		$map_file = fastq2LASTpe() if($last_dir);
		$map_file = fastq2GEMpe() if($gem_dir);
	}else{
		$map_file = fastq2SOAPse() if($soap2_exe);
		$map_file = fastq2BOWTIEse() if($bowtie2_exe);
		$map_file = fastq2BWAse() if($bwa_exe);
		$map_file = fastq2LASTse() if($last_dir);
		$map_file = fastq2GEMse() if($gem_dir);
	}
	my $end_time = time;
	my $time_taken = $end_time - $start_time;	
	print "Finished mapping and converting reads in $time_taken\n";
	sortsam($map_file);
	undef %chrSizes;
	my $total_mreads = 0;
	foreach my $chr(keys %countReads){
		print $chr, "\t", $countReads{$chr}, "\n";
		$total_mreads += $countReads{$chr};
	}
	print "Total sequences\t$total_reads\n";
	print "Total bases\t$total_bases\n";
	print "Total mapped bases\t$total_mbases\n";
	print "Total wrong guesses\t", $missed_guesses, "\n";
	undef %countReads;
}

sub processhit{
	my $array = shift;
	my @f = @{$array};
	my ($id, $orig_seq, $guess) = split /\|/, $f[0];
	$f[0] = $id;
	$f[9] = $orig_seq;
	my @qual = split "", $f[10];
	if($f[2] =~ s/_Watson//){
		if($f[1] & 16){
			return if($guess ne "R");
			# make everything on Watson maps to forward
			$f[9] = revComp($orig_seq);
			$f[10] = join("", reverse(@qual));
			#$f[0] = $f[0].":R";
		}else{
			return if($guess ne "F");
			#$f[0] = $f[0].":F";
		}
		$f[1] = 0;
	}elsif($f[2] =~ s/_Crick//){
		if($f[1] & 16){
			return if($guess ne "F");
			#$f[0] = $f[0].":R";
			# reverse on Crick is forward on Watson, make it reverse on Watson
			$f[9] = revComp($orig_seq);
			$f[10] = join("", reverse(@qual));
		}else{
			return if($guess ne "R");
			#$f[0] = $f[0].":F";
		}
		$f[1] = 16;
	}
	###-----------------BEGIN deal with CIGAR---------------------###
	my @CIGARS;
		# I - insertion to the reference
		# D - deletion to the reference
		# N - skipped region from the reference
		# S - soft clip on the read
		# H - hard clip on the read
		# P - padding - silent deletion from the padded reference sequence
	my @values = split /(\d+)/, $f[5];
	my $hard_clip_beg = 0;
	my $hard_clip_end = 0;
	for(my $i = 1; $i<scalar(@values)-1; $i=$i+2){
		if($values[$i+1] eq "H"){ #hide hard clippings
			$hard_clip_beg = $values[$i] if($i == 1);
			$hard_clip_end = $values[$i] if($i > 1);
			next;
		}
		push(@CIGARS, $values[$i].$values[$i+1]);
	}
	# match orig-seq to CIGARS
	my $clipped_s = substr($f[9], $hard_clip_beg, length($f[9])-$hard_clip_beg-$hard_clip_end);
	$f[5] = join("", @CIGARS);
	$f[9] = $clipped_s;
	###-----------------END deal with CIGAR---------------------###

	$total_mbases+=length($f[9]);
	$countReads{$f[2]}++;
	return @f[0...10];
}

sub sortsam{
	my $map_file = shift;
	#identify unique reads and save them in two files based on the templates.
	my $combined_sorted_map_file = $fqName.".combined.sorted";
	my $cmd = "sort -k 1,1 -T $TMP_DIR < $map_file > $combined_sorted_map_file";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	unlink($map_file);
	
	open(SAM_IN, "$combined_sorted_map_file") || die("Error in opening $combined_sorted_map_file.");

	my %fileHandle;
	foreach my $chr(keys %chrSizes){
		my $fname = $fqName . "." . $chr . ".sam";
		open ( $fileHandle{$chr}->{"h"} , ">$fname") || die("Error writing to file $fname\n");
		$fileHandle{$chr}->{"n"} = $fname;
	}
	my $last_line = 0;
	my @last_fields;
	my $last_cnt = 0;
	while(my $line =  <SAM_IN>){
		#next if($line =~ m/XS:i/ && $bowtie2_exe);
		my @fields = split(/\t/, $line);
		next if(scalar(@fields) < 5); # skip the header lines
		#next if($fields[4] < 5); # skip if MAPQ is less than 5
		if($fields[5]){ # CIGAR
			next if($fields[5] eq "*");
			my @cigar_str = split(/(\d+)/, $fields[5]);
			my $identity = 0;
			for(my $i = 0; $i < scalar(@cigar_str); $i++){
				if($cigar_str[$i] eq "M"){
					$identity+= $cigar_str[$i-1];
				}
			}
			next if($identity/length($fields[9]) < $min_identity); # skip read if less than $min_identity
		}
		if(!$last_line){
			$last_line = $line;
			@last_fields = @fields;
			$last_cnt = 0;
			next;
		}
		if($fields[0] eq $last_fields[0]){
			#my @tmp = split ":", $fields[11];
			#my $score = pop(@tmp);
			#@tmp = split ":", $last_fields[11];
			#my $last_score = pop(@tmp);
			#undef(@tmp);
			my $score = $fields[4];
			my $last_score = $last_fields[4];
			if($score > $last_score){
				#2nd line is a better hit
				$last_line = $line;
				@last_fields = @fields;
				$last_cnt = 0;
			}elsif($score == $last_score){
				#two equivalent good hits, increment last_cnt
				$last_cnt++;
			}else{
				#1st line is a better hit, do nothing
			}
		}else{	
			if($last_cnt ne 0){ # not a unique best hit
				undef $last_line;
				undef @last_fields;
				next;
			}
			if(my @hit = processhit(\@last_fields)){
				print {$fileHandle{$hit[2]}->{"h"}} join("\t", @hit), "\n";
			}else{
				print "Wrong guess for ", $last_fields[0], "\n";
				$missed_guesses++;
			}
			$last_line = $line;
			@last_fields = @fields;
		}
	}
	#print the last line
	if($last_line and $last_cnt eq 0){
		if(my @hit = processhit(\@last_fields)){
			print {$fileHandle{$hit[2]}->{"h"}} join("\t", @hit), "\n";
		}else{
			print "Wrong guess for ", $last_fields[0], "\n";
			$missed_guesses++;
		}
	}
	close(SAM_IN);
	unlink($combined_sorted_map_file);
	foreach my $chr (keys %fileHandle){
		close($fileHandle{$chr}->{"h"});
		my $unsorted = $fileHandle{$chr}->{"n"};
		$cmd = "sort -T $TMP_DIR -k4,4n $unsorted > $fqName.$chr.sorted.sam";
		system($cmd) == 0 or die "system problem (exit $?): $!\n";
		unlink($unsorted);
	}
}

###----------------- SOAP mapper-------------------###
sub fastq2SOAPse(){
	my $soap_fwd_map_file = $fqName.".fwd.soap.out";
	my $cmd = "$soap2_exe -r 0 -v $maxMismatches -p $cpu -D $template_fwd -a $encodedFqName1 -o $soap_fwd_map_file";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	print "$cmd\n";
        my $map_file = $fqName . ".soap.sam";
        soap2sam($soap_fwd_map_file, 0, $map_file);
	unlink($soap_fwd_map_file);
	if($template_rev){
		my $soap_rev_map_file = $fqName.".rev.soap.out";
		$cmd = "$soap2_exe -r 0 -v $maxMismatches -p $cpu -D $template_rev -a $encodedFqName1 -o $soap_rev_map_file";
		system($cmd) == 0 or die "system problem (exit $?): $!\n";
		print "$cmd\n";
		soap2sam($soap_rev_map_file, 0, $map_file);
		unlink($soap_rev_map_file);
	}
	return ($map_file);
}

sub fastq2SOAPpe(){
	my $soap_fwd_map_PE_file = $fqName.".fwd.soap.PE.out";
	my $soap_fwd_map_SE_file = $fqName.".fwd.soap.SE.out";
	my $cmd = "$soap2_exe -r 0 -v $maxMismatches -p $cpu -D $template_fwd -a $encodedFqName1 -b $encodedFqName2 -o $soap_fwd_map_PE_file -2 $soap_fwd_map_SE_file -m 1 -x 1000";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	print "$cmd\n";
        my $map_file = $fqName . ".soap.PE.sam";
        soap2sam($soap_fwd_map_PE_file, 0, $map_file);
        soap2sam($soap_fwd_map_SE_file, 0, $map_file);
        unlink($soap_fwd_map_SE_file);
        unlink($soap_fwd_map_PE_file);

	if($template_rev){
		my $soap_rev_map_PE_file = $fqName.".rev.soap.PE.out";
		my $soap_rev_map_SE_file = $fqName.".rev.soap.SE.out";
		$cmd = "$soap2_exe -r 0 -v $maxMismatches -p $cpu -D $template_rev -a $encodedFqName1 -b $encodedFqName2 -o $soap_rev_map_PE_file -2 $soap_rev_map_SE_file -m 1 -x 1000";
		system($cmd) == 0 or die "system problem (exit $?): $!\n";
		print "$cmd\n";
		soap2sam($soap_rev_map_PE_file, 0, $map_file);
		soap2sam($soap_rev_map_SE_file, 0, $map_file);
		unlink($soap_rev_map_SE_file);
		unlink($soap_rev_map_PE_file);
	}

	return ($map_file);
}
###----------------- END SOAP mapper-------------------###

###----------------- BOWTIE2 mapper-------------------###
sub fastq2BOWTIEse(){
	# Set the correct mapping parameters
	my $options = "--phred$qual_base --very-sensitive --score-min $score_min -p $cpu --no-unal --no-head";
	my $map_file = $fqName.".bowtie.sam";	
	my $cmd = "$bowtie2_exe $options -x $template_fwd -U $encodedFqName1 | grep -v XS:i: > $map_file";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	print "$cmd\n";

	return ($map_file) if(!$template_rev);
	$cmd = "$bowtie2_exe $options -x $template_rev -U $encodedFqName1 | grep -v XS:i: >> $map_file"; 
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	print "$cmd\n";
	return ($map_file);
}
sub fastq2BOWTIEpe(){
	# Set the correct mapping parameters
	my $options = "--phred$qual_base --very-sensitive --score-min $score_min -p $cpu -I 0 -X 1000 --no-discordant --no-unal --no-head";
	my $map_file = $fqName.".bowtie.PE.sam";
	# note: XS:i: defined as multimapping is specific to bowtie only.
	my $cmd = "$bowtie2_exe $options -x $template_fwd -1 $encodedFqName1 -2 $encodedFqName2 | grep -v XS:i: > $map_file";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	print "$cmd\n";

	return ($map_file) if(!$template_rev);
	$cmd = "$bowtie2_exe $options -x $template_rev -1 $encodedFqName1 -2 $encodedFqName2 | grep -v XS:i: >> $map_file";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	print "$cmd\n";
	
	return ($map_file);
}
###----------------- END BOWTIE2 mapper-------------------###

###----------------- BWA mem mapper-------------------###
sub fastq2BWAse(){
        # Set the correct mapping parameters
        my $options = "mem -t $cpu -B2 -c 1000";
        my $map_file = $fqName.".bwa.sam";

        my $cmd = "$bwa_exe $options $template_fwd $encodedFqName1 | awk '{if(\$5 > 5) print \$0}' > $map_file";
        system($cmd) == 0 or die "system problem (exit $?): $!\n";
        print "$cmd\n";

	return ($map_file) if(!$template_rev);
        $cmd = "$bwa_exe $options $template_rev $encodedFqName1 | awk '{if(\$5>5) print \$0}' >> $map_file";
        system($cmd) == 0 or die "system problem (exit $?): $!\n";
        print "$cmd\n";
	return ($map_file);
}
sub fastq2BWApe(){
        # Set the correct mapping parameters
        my $options = "mem -t $cpu -B2 -c 1000 -a";
        my $map_file = $fqName.".bwa.PE.sam";

        my $cmd = "$bwa_exe $options $template_fwd $encodedFqName1 $encodedFqName2 | awk '{if(\$5>5) print \$0}' > $map_file";
        system($cmd) == 0 or die "system problem (exit $?): $!\n";
        print "$cmd\n";
	
	return ($map_file) if(!$template_rev);
        $cmd = "$bwa_exe $options $template_rev $encodedFqName1 $encodedFqName2 | awk '{if(\$5>5)print \$0}' >> $map_file";
        system($cmd) == 0 or die "system problem (exit $?): $!\n";
        print "$cmd\n";
        return ($map_file);
}
###----------------- END BWA mem mapper-------------------###

###----------------- LAST mapper-------------------###
sub fastq2LASTse(){
        my $last_map_probs = "$last_dir/scripts/last-map-probs.py -s150";
        my $lastal = "$last_dir/src/lastal -Q1 -e120";
        my $maf_convert = "$last_dir/scripts/maf-convert.py sam";
        my $map_file = $fqName.".last.SE.sam";
        my $cmd = "$lastal $template_fwd $encodedFqName1 | $last_map_probs - | $maf_convert > $map_file";
        system($cmd) == 0 or die "system problem (exit $?): $!\n";
        print "$cmd\n";
        return ($map_file) if(!$template_rev);
        $cmd = "$lastal $template_rev $encodedFqName1 | $last_map_probs - | $maf_convert >> $map_file";
        system($cmd) == 0 or die "system problem (exit $?): $!\n";
        print "$cmd\n";
        return ($map_file);
}
sub fastq2LASTpe(){
	# NOTE: last_pair_probs.py is not used to save hard disk requirements (pipe results to last_map_probs!)
	my $last_map_probs = "$last_dir/scripts/last-map-probs.py -s150";
	my $lastal = "$last_dir/src/lastal -Q1 -e120";
	my $maf_convert = "$last_dir/scripts/maf-convert.py sam";
	my $map_file = $fqName.".last.PE.sam";
	my $cmd = "$lastal $template_fwd $encodedFqName1 $encodedFqName2 | $last_map_probs - | $maf_convert > $map_file";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	print "$cmd\n";
	return ($map_file) if(!$template_rev);
	$cmd = "$lastal $template_rev $encodedFqName1 $encodedFqName2 | $last_map_probs - | $maf_convert >> $map_file";
	system($cmd) == 0 or die "system problem (exit $?): $!\n";
	print "$cmd\n";
	return ($map_file);
}
###----------------- END LAST mapper-------------------###

###----------------- GEM mapper-------------------###
sub fastq2GEMse(){
        # Set the correct mapping parameters
	my $gem_exe="$gem_dir/gem-mapper";
	my $gem_2_sam="$gem_dir/gem-2-sam";
	my $offset = "offset-33";
	$offset = "offset-64" if($qual_base eq 64);
        my $options = "-i $encodedFqName1 -q $offset -T $cpu $gem_allow";
        my $map_file = $fqName.".gem.sam";

        my $cmd = "$gem_exe -I $template_fwd $options | $gem_2_sam -I $template_fwd -q $offset -c | awk '{if(\$3 != \"*\") print \$0}' > $map_file";
        system($cmd) == 0 or die "system problem (exit $?): $!\n";
        print "$cmd\n";

	return ($map_file) if(!$template_rev);
        $cmd = "$gem_exe -I $template_rev $options | $gem_2_sam -I $template_rev -q $offset -c | awk '{if(\$3 != \"*\") print \$0}' >> $map_file";
        system($cmd) == 0 or die "system problem (exit $?): $!\n";
        print "$cmd\n";
	return ($map_file);
}
sub fastq2GEMpe(){
        # Set the correct mapping parameters
	my $gem_exe="$gem_dir/gem-mapper";
	my $gem_2_sam="$gem_dir/gem-2-sam";
	my $offset = "offset-33";
	$offset = "offset-64" if($qual_base eq 64);
	# Pair end mapping with GEM is currently not possible.
        my $options = "-q $offset -T $cpu $gem_allow";
        my $map_file = $fqName.".gem.PE.sam";
	print "Paired-end mapping with GEM is currently unavailable\n";

        my $cmd = "$gem_exe -I $template_fwd -i $encodedFqName1 $options | $gem_2_sam -I $template_fwd -q $offset -c | awk '{if(\$3 != \"*\") print \$0}' > $map_file";
        system($cmd) == 0 or die "system problem (exit $?): $!\n";
        print "$cmd\n";
        $cmd = "$gem_exe -I $template_fwd -i $encodedFqName2 $options | $gem_2_sam -I $template_fwd -q $offset -c | awk '{if(\$3 != \"*\") print \$0}' >> $map_file";
        system($cmd) == 0 or die "system problem (exit $?): $!\n";
        print "$cmd\n";

	return ($map_file) if(!$template_rev);

        $cmd = "$gem_exe -I $template_rev -i $encodedFqName1 $options | $gem_2_sam -I $template_rev -q $offset -c | awk '{if(\$3 != \"*\") print \$0}' >> $map_file";
        system($cmd) == 0 or die "system problem (exit $?): $!\n";
        print "$cmd\n";
        $cmd = "$gem_exe -I $template_rev -i $encodedFqName2 $options | $gem_2_sam -I $template_rev -q $offset -c | awk '{if(\$3 != \"*\") print \$0}' >> $map_file";
        system($cmd) == 0 or die "system problem (exit $?): $!\n";
        print "$cmd\n";

	return ($map_file);
}
###----------------- END GEM mapper-------------------###

sub encodeFastq{
	my $fileName = shift;
	my $outFileName = shift;
	my $fivep = $trim5;
	my $threep = $trim3;
	if($fileName =~ /\.gz$/){
		open(FQ, "gunzip -c $fileName |") || die("Can't open pipe to $fileName!");
	}else{
		open(FQ, "$fileName") || die ("Error in opening file $fileName!");
	}
	open(FQ_OUT, ">$outFileName") || die ("Error in opening file $outFileName!");
	while(my $line1 = <FQ>){
		chomp($line1);
		my @tmp = split /[\s+\t]/, $line1;
		# first 6 fields separated by : are cluster ID and by _ are UMI
		$line1 = $tmp[0];
		my $line2 = <FQ>; 
		$line2 =~ tr/\./N/;
		chomp($line2);
		my $line3 = <FQ>;
		chomp($line3);
		my $line4 = <FQ>;
		chomp($line4);
		next if(!$line1 || !$line2 || !$line3 || !$line4);
		if($threep || $fivep){
			my $start = $fivep;
			my $total = length($line2) - $threep - $fivep;
			next if($total < 20);
			$line2 = substr($line2, $start, $total);
			$line4 = substr($line4, $start, $total);
		}
		if($qualtrim){
		#http://wiki.bioinformatics.ucdavis.edu/index.php/TrimBWAstyle.pl
			my @quals = split "", $line4;
			my $pos = scalar(@quals);
			my $maxPos = $pos;
			my $sum_bad_qual = 0;
			my $maxSum = 0;
			while( $pos > 0 && $sum_bad_qual >= 0){
				$sum_bad_qual += $qualtrim - (ord($quals[$pos-1])-$qual_base);
				if($sum_bad_qual > $maxSum){
					$maxSum = $sum_bad_qual;
					$maxPos = $pos;
				}
				$pos--;
			}
			if($pos==0 || $maxPos < 20) { $line2 = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"; 
						      $line4 = "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"; }
			else{
				$line2 = substr($line2,0,$maxPos);
				$line4 = substr($line4,0,$maxPos);
			}
		}
		my $guess = guess_strand($line2);
		my $prefix = $line1;
		die("Read names are too long!\n$line1\n") if(length($prefix) > 250 - 7);
		my $s = $line2;
		if($guess eq "R"){
			$line2 =~ s/G/a/g;
		}else{
			$line2 =~ s/C/t/g;
		}
		$line1 = $prefix . "|" . $s . "|" . $guess;
		if(length($line1) > 250){
			#guess takes 2 characters
			#number takes 5 characters
			my $max_length = 250 - length($prefix) - 7;
			# try not to cut into pieces that are too small
			my $j = $max_length;
			while(length($s) % $j < 50 && $j > 50){
				$j--;
			}
			$max_length = $j;
			my $leftoverOrig = $s;
			my $leftoverNew = $line2;
			my $leftoverQual = $line4;
			my $i = 0;
			my $cnt = 0;
			while($i < length($leftoverOrig)){
				my $l = length($leftoverOrig) - $i;
				$l = $max_length if(length($leftoverOrig) - $i >= $max_length);
				$s = substr($leftoverOrig, $i, $l);
				$line2 = substr($leftoverNew, $i, $l);
				$line4 = substr($leftoverQual, $i, $l);
				$i+=$l;
				$cnt++;
				$line1 = $prefix . "_" . $cnt . "|" . $s . "|" . $guess;
				$total_bases+=length($line2);
				print FQ_OUT "$line1\n$line2\n$line3\n$line4\n";
			}
		}else{
			$total_bases+=length($line2);
			$total_reads++;	
			print FQ_OUT "$line1\n$line2\n$line3\n$line4\n";
		}
	}
	close(FQ);
	close(FQ_OUT);
}

sub revComp{
	my $seq = shift;
	my $rcSeq='';
	for(my $i=0; $i<length($seq); $i++){
		$rcSeq = $rcTable{uc(substr($seq,$i,1))} . $rcSeq;
	}
	return $rcSeq;
}

sub guess_strand{
	my $seq = shift;
	my %baseCounts;
	$baseCounts{'A'}=0.001;
	$baseCounts{'T'}=0.001;
	$baseCounts{'G'}=0.001;
	$baseCounts{'C'}=0.001;
	while(my $base = chop($seq)){
		$baseCounts{$base}++;
	}
	if($baseCounts{'T'}/$baseCounts{'C'} > $baseCounts{'A'}/$baseCounts{'G'}) {
		return "F";
	}else{
		return "R";
	}
}

sub getChromSizes{
	return 0 if($template_idx !~ m/fai/);
	open(GENOME_INDEX, "$template_idx") || return 0;
	while(my $line = <GENOME_INDEX>){
		chomp($line);
		my @f = split "\t", $line;
		my $cur_chr = $f[0];
		$cur_chr =~ s/_Watson//;
		$cur_chr =~ s/_Crick//;
		$f[0] = $cur_chr;
		$chrSizes{$cur_chr} = $f[1];
	}
	close(GENOME_INDEX);
	return 1;
}

sub printUsage(){
	print "Usage: bisReadMapper.pl [options] > log\n";
	print "Required options:\n";
	print "   -r   read(s) fq1[,fq2]\n";
	print "   -W   path to Watson converted index OR single combined W/C index\n";
	print "   -C   path to Crick converted index\n";
	print "   -g   path to indexed reference genome fasta file (.fai) \n";
	print "   -a   path to aligner, SOAP2 or Bowtie2 or path to LAST directory or path to GEM binary directory only\n";
	print "   -e   encoded [1/0]\n";
	print "        If yes, then reads will not be encoded nor trimmed\n";
	print "   -p   number of processors for mapping [int], default 2\n";
	print "   -n   name of sample\n";
	print "   -3   3-prime trimming\n";
	print "   -b   ASCII base quality offset\n";
	print "   -5   5-prime trimming\n";
	print "   -q   quality value to perform quality trimming at\n";
	print "   -T   temporary directory for unix sort\n";
	print "   -s   --min-score function for bowtie2 (-L,-0.6,-0.6)\n";
}

###-------- soap2sam.pl from heng li--------------###

sub mating {
  my ($s1, $s2) = @_;
  my $isize = 0;
  if ($s1->[2] ne '*' && $s1->[2] eq $s2->[2]) { # then calculate $isize
	my $x1 = ($s1->[1] & 0x10)? $s1->[3] + length($s1->[9]) : $s1->[3];
	my $x2 = ($s2->[1] & 0x10)? $s2->[3] + length($s2->[9]) : $s2->[3];
	$isize = $x2 - $x1;
  }
  # update mate coordinate
  if ($s2->[2] ne '*') {
	@$s1[6..8] = (($s2->[2] eq $s1->[2])? "=" : $s2->[2], $s2->[3], $isize);
	$s1->[1] |= 0x20 if ($s2->[1] & 0x10);
  } else {
	$s1->[1] |= 0x8;
  }
  if ($s1->[2] ne '*') {
	@$s2[6..8] = (($s1->[2] eq $s2->[2])? "=" : $s1->[2], $s1->[3], -$isize);
	$s2->[1] |= 0x20 if ($s1->[1] & 0x10);
  } else {
	$s2->[1] |= 0x8;
  }
}

sub soap2sam {
  my ($in_file, $is_paired, $out_file) = @_;
  # core loop
  my @s1 = ();
  my @s2 = ();
  my ($s_last, $s_curr) = (\@s1, \@s2);
  open(INFILE, "$in_file") || die("Error opening $in_file\n");
  open(OUTFILE, ">>$out_file") || die("Error writing to $out_file\n");
  while (<INFILE>) {
	s/[\177-\377]|[\000-\010]|[\012-\040]//g;
	next if (&soap2sam_aux($_, $s_curr, $is_paired) < 0);
	if (@$s_last != 0 && $s_last->[0] eq $s_curr->[0]) {
	  &mating($s_last, $s_curr);
	  print OUTFILE join("\t", @$s_last), "\n";
	  print OUTFILE join("\t", @$s_curr), "\n";
	  @$s_last = (); @$s_curr = ();
	} else {
	  print OUTFILE join("\t", @$s_last), "\n" if (@$s_last != 0);
	  my $s = $s_last; $s_last = $s_curr; $s_curr = $s;
	}
  }
  print OUTFILE join("\t", @$s_last), "\n" if (@$s_last != 0);
  close(INFILE);
  close(OUTFILE);
}

sub soap2sam_aux {
  my ($line, $s, $is_paired) = @_;
  chomp($line);
  my @t = split(/\s+/, $line);
  return -1 if (@t < 9 || $line =~ /^\s/ || !$t[0]);
  @$s = ();
  # fix SOAP-2.1.x bugs
  @t = @t[0..2,4..$#t] unless ($t[3] =~ /^\d+$/);
  # read name
  $s->[0] = $t[0];
  $s->[0] =~ s/\/[12]$//g;
  # initial flag (will be updated later)
  $s->[1] = 0;
  $s->[1] |= 1 | 1<<($t[4] eq 'a'? 6 : 7);
  $s->[1] |= 2 if ($is_paired);
  # read & quality
  $s->[9] = $t[1];
  $s->[10] = (length($t[2]) > length($t[1]))? substr($t[2], 0, length($t[1])) : $t[2];
  # cigar
  $s->[5] = length($s->[9]) . "M";
  # coor
  $s->[2] = $t[7]; $s->[3] = $t[8];
  $s->[1] |= 0x10 if ($t[6] eq '-');
  # mapQ
  $s->[4] = $t[3] == 1? 30 : 0;
  # mate coordinate
  $s->[6] = '*'; $s->[7] = $s->[8] = 0;
  # aux
  push(@$s, "NM:i:$t[9]");
  my $md = '';
  if ($t[9]) {
	my @x;
	for (10 .. $#t) {
	  push(@x, sprintf("%.3d,$1", $2)) if ($t[$_] =~ /^([ACGT])->(\d+)/i);
	}
	@x = sort(@x);
	my $a = 0;
	for (@x) {
	  my ($y, $z) = split(",");
	  $md .= (int($y)-$a) . $z;
	  $a += $y - $a + 1;
	}
	$md .= length($t[1]) - $a;
  } else {
	$md = length($t[1]);
  }
  push(@$s, "MD:Z:$md");
  return 0;
}
