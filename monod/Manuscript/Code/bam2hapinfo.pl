#!/usr/bin/perl -w

use strict;
use warnings;
use Cwd;
use Getopt::Long;
use Pod::Usage;
use IO::Handle;
use IO::Socket::INET;

my $parent_dir = getcwd();
my $start_run = time();
warn "# Run started at: $start_run\n";
my $command_line = join (" ",@ARGV);
my $samtools = "samtools";
my ($merged_bam_file,$target_list_file,$library,$cpg_position_file,
           $chrsizeFile,$Aligner,$phred,$server,$queue,$help,$pbs,$submit,$nodes,$ppn,
           $walltime,$multicore)=&process_command_line();

if($pbs){
	&pbsprocess(\$merged_bam_file,\$target_list_file,\$library,\$cpg_position_file,
	\$chrsizeFile,\$Aligner,\$phred,\$queue,\$ppn,\$nodes,\$walltime,\$parent_dir);
	exit;
}


my $target_flanking_len = 80;
my %targetTable;                     # stores target regions          
my %hapInfoTable;                    # stores haplotype
my %binnedtargetTable;               # stores genomic bins, cpgs and genomic regions
my %rcTable;                         # stores rctable
my $verbose = 0;
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
my $reads_parsed ;                    # counts the reads number in the bam files
my $reads_after_filtering;            # counts the reads number used for hapinfo
open F,"$chrsizeFile" || die "can not open $chrsizeFile,please make sure this file is existed!";
my %chrSizes;
while(<F>){
chomp;
my($chr,$len)=split/\s+/;
$chrSizes{$chr}=$len;
}
# Step 1. Assign target region(targetTable) to genomic bins (binnedtargetTable)
open(INFILE, "$target_list_file")||die("Error in opening $target_list_file\n");
my $bin_size = 100000; # warn: keep your interest regions < $bin_size
while(my $line = <INFILE>){
	next if $line !~/^chr/;
	chomp($line);
        next if $line=~/^\s+$/;
	my ($chr,$target_start,$target_end,$id) = split(/\s+/, $line);
	my $len=$target_end-$target_start;
    warn "# $chr:$target_start-$target_end is longer than $bin_size bp, please contact Shihcheng.Guo\@gmail.com\n\n" if $len > $bin_size;
	my $index = int($target_start/$bin_size);
	# print "index :$index\n";                        #  Debug(20170201) 
	my $target_id = "$chr:$target_start-$target_end";
	push(@{$binnedtargetTable{$chr}->{$index}},$target_id);
	push(@{$binnedtargetTable{$chr}->{$index-1}},$target_id);
	push(@{$binnedtargetTable{$chr}->{$index+1}},$target_id);
}
close(INFILE);

# Step 2. Assign CpG loci to genomic bin (binnedtargetTable) and target region( targetTable)
open(INFILE, "$cpg_position_file")|| die "Error in opening $cpg_position_file, Please make sure this file existed!\n";
while(my $line = <INFILE>){
	chop($line);
	my ($chr, $pos) = split(/\s+/, $line);
	my $index = int($pos/$bin_size);
	next if(!$binnedtargetTable{$chr}->{$index});	
	# print "index :$index\n";                     # Debug(20170201) 
	foreach my $target_id (@{$binnedtargetTable{$chr}->{$index}}){
		my ($chr,$target_start,$target_end) = split(/[:\-]/, $target_id);
		if($target_start <= $pos && $pos <=$target_end){
		push(@{$targetTable{$target_id}->{"CpG_positions"}},$pos);	
		}
	}
}
close(INFILE);

# Step 3. loop interest regions and parse CpG methylation status with samtools view. 
if($Aligner eq 'bismark'){
foreach my $target_id (sort keys(%targetTable)){
	my ($chr,$target_start,$target_end) = split(/[:\-]/, $target_id);
	my $region_start = $target_start-$target_flanking_len;
	my $region_end = $target_end+$target_flanking_len;
	my $cmd = "$samtools view -q 10 $merged_bam_file $chr:$region_start-$region_end";
	my $ans = `$cmd`;
	my @lines = split(/\n/, $ans);
	next if(scalar(@lines)<1 || !$targetTable{$target_id}->{"CpG_positions"});
	$reads_parsed += scalar(@lines);
	$reads_after_filtering += scalar(@lines);
	my @CpG_positions = sort {$a <=> $b} @{$targetTable{$target_id}->{"CpG_positions"}};
	my %readHapInfo;	

	#print "##$cmd\n";
	#Going through the reads one at a time
	my %read_start_pos_table;
	foreach my $line (@lines){
		my @fields = split(/\t/, $line);
		my $read_strand = $fields[1];		
		my $seq = $fields[9];		
		my $qual_string = $fields[10];
		my $read_start = $fields[3];
		my $CIGAR = $fields[5];
		my $read_length = length($seq);
		next if($CIGAR =~ /[ID]/);
		
		if($CIGAR=~ /S/){
			my ($clip_len, @others) = split(/S/, $CIGAR);
			next if(length($clip_len)>2);
			$seq=substr($seq,$clip_len,$read_length-$clip_len);		
			$qual_string=substr($qual_string,$clip_len,$read_length-$clip_len);	
			$read_length -= $clip_len;			
		}elsif($seq =~ /^CG/){
			$seq=substr($seq,2,$read_length-2);		
			$qual_string=substr($qual_string,2,$read_length-2);	
			$read_start += 2;
			$read_length -= 2;
		}
		my @read_fields = split(/[\:\#\_]+/, $fields[0]);
		my $N_read_fields = scalar(@read_fields);
		for(my $i=$N_read_fields; $i<=4; $i++){
			push(@read_fields,"NA");
		}
		my $read_id = scalar(@read_fields)>7 ? join("-", @read_fields[0..6]) : join("-", @read_fields[0..4]);	
		#print "#$read_id\t$CIGAR\t$read_start\t$read_strand\t$seq\t$qual_string\n";		
		foreach my $CpG_position(@CpG_positions){
			my $offset = $CpG_position-$read_start+1;
			next if($offset<0 || $offset >= length($seq));	 # offset should not or cannot <0 or else you should change the code or data carefully.
			# Debug: print "$CpG_position\t$read_start\t$offset\n";
            # the situation sometimes should be change dependent on Flag of different alignmentor
			if($read_strand <16 || $read_strand ==99 ||$read_strand ==147)
            {
                # positive chain
				my $qual_score = ord(substr($qual_string,$offset,1))-$phred;
				next if($readHapInfo{$read_id}->{$CpG_position}->{"qual"} && $readHapInfo{$read_id}->{$CpG_position}->{"qual"} > $qual_score); # UMIs choose the best one
				$readHapInfo{$read_id}->{$CpG_position}->{"base"}=substr($seq,$offset,1);
				$readHapInfo{$read_id}->{$CpG_position}->{"qual"}=$qual_score;
			}else{
				# Negative chain
				my $qual_score = ord(substr($qual_string,$offset+1,1))-$phred;
				next if($readHapInfo{$read_id}->{$CpG_position}->{"qual"} && $readHapInfo{$read_id}->{$CpG_position}->{"qual"} > $qual_score); # UMIs choose the best one
				$readHapInfo{$read_id}->{$CpG_position}->{"base"}=$rcTable{substr($seq,$offset+1,1)};			
				$readHapInfo{$read_id}->{$CpG_position}->{"qual"}=$qual_score;	

			}	
			#print "$offset\t", $readHapInfo{$read_id}->{$CpG_position}, "\n";
		}	
		$read_start_pos_table{$read_id} = $read_start if(!$read_start_pos_table{$read_id} || $read_start_pos_table{$read_id}<$read_start);
	}		
	
	# Use read start and read end as the UMI to collapse multiple clonal reads.
	# How to define UMI will influence the result compared with classic methylation level(bismark)
	my @read_list = keys(%readHapInfo);		
	my %unique_read_base_info;
	foreach my $read_id (@read_list){
	my $UMI;
	if($library=="RRBS"){
	$UMI = $read_id; # read_id as UMI 
	}else{                        
	$UMI = $read_start_pos_table{$read_id};  # start postion as UMI
	}
	foreach my $CpG_position (keys(%{$readHapInfo{$read_id}})){
		next if(!$readHapInfo{$read_id}->{$CpG_position}->{"base"});
		#print "$CpG_position\n";
		#print $readHapInfo{$read_id}->{$CpG_position}->{"base"}, "\n";
		#print $readHapInfo{$read_id}->{$CpG_position}->{"qual"}, "\n";
		$unique_read_base_info{$UMI}->{$CpG_position}->{$readHapInfo{$read_id}->{$CpG_position}->{"base"}}+=$readHapInfo{$read_id}->{$CpG_position}->{"qual"};
		}
	}
	
	#Derive the consensus haplotype string based on multiple clonal reads.
	foreach my $UMI (keys(%unique_read_base_info)){
		my $hap_string="";
		foreach my $CpG_position(@CpG_positions){
			my $base = "N";
			my $best_qual = 0;
			if($unique_read_base_info{$UMI}->{$CpG_position}){
				foreach my $base_call (keys(%{$unique_read_base_info{$UMI}->{$CpG_position}})){
					next if($unique_read_base_info{$UMI}->{$CpG_position}->{$base_call} <= $best_qual);
					$best_qual = $unique_read_base_info{$UMI}->{$CpG_position}->{$base_call};
					$base = $base_call;
				}
			}					
			$hap_string=$hap_string.$base;
		}	
		$hapInfoTable{$target_id}->{"hap_counts"}->{$hap_string}++ if($hap_string);
	}

	#report haplotype strings only on the positions with valid calls.
	foreach my $hap_string(keys %{$hapInfoTable{$target_id}->{"hap_counts"}}){
		my @valid_positions;
		my $valid_hap="";
		for(my $i=0; $i<length($hap_string); $i++){
			my $allele = substr($hap_string,$i,1);
			next if($allele eq "N");
			$valid_hap = $valid_hap.$allele;
			push(@valid_positions,$CpG_positions[$i]);
		}
		next if(scalar(@valid_positions)<1);
		next if($valid_hap =~ /[AG]/);
		print $target_id,"\t$valid_hap\t", $hapInfoTable{$target_id}->{"hap_counts"}->{$hap_string},"\t",join(",",@valid_positions), "\n";
	}
}
}else{

foreach my $target_id (sort keys(%targetTable)){
	my ($chr,$target_start,$target_end) = split(/[:\-]/, $target_id);
	my $region_start = $target_start-$target_flanking_len;
	my $region_end = $target_end+$target_flanking_len;
	my $cmd = "$samtools view -q 10 $merged_bam_file $chr:$region_start-$region_end";
	my $ans = `$cmd`;
	my @lines = split(/\n/, $ans);
	next if(scalar(@lines)<1 || !$targetTable{$target_id}->{"CpG_positions"});
	$reads_parsed += scalar(@lines);
	$reads_after_filtering += scalar(@lines);
	my @CpG_positions = sort {$a <=> $b} @{$targetTable{$target_id}->{"CpG_positions"}};   # @CpG_positions contains all the CpGs within target_id
	my %readHapInfo;	

	#print "##$cmd\n";
	#Going through the reads one at a time
	my %read_start_pos_table;
	foreach my $line (@lines){
		my @fields = split(/\t/, $line);
		my $read_strand = $fields[1];		
		my $seq = $fields[9];		
		my $qual_string = $fields[10];
		my $read_start = $fields[3];
		my $CIGAR = $fields[5];
		my $read_length = length($seq);
		next if($CIGAR =~ /[ID]/);
        
		if($CIGAR=~ /S/){
			my ($clip_len, @others) = split(/S/, $CIGAR);
			next if(length($clip_len)>2);
			$seq=substr($seq,$clip_len,$read_length-$clip_len);		
			$qual_string=substr($qual_string,$clip_len,$read_length-$clip_len);	
			$read_length -= $clip_len;			
		}elsif($seq =~ /^CG/){
			$seq=substr($seq,2,$read_length-2);		
			$qual_string=substr($qual_string,2,$read_length-2);	
			$read_start += 2;
			$read_length -= 2;
		}
		my @read_fields = split(/[\:\#\_]+/, $fields[0]);
		my $N_read_fields = scalar(@read_fields);
		for(my $i=$N_read_fields; $i<=4; $i++){
			push(@read_fields,"NA");
		}
		my $read_id = scalar(@read_fields)>7 ? join("-", @read_fields[0..6]) : join("-", @read_fields[0..4]);	
		#print "#$read_id\t$CIGAR\t$read_start\t$read_strand\t$seq\t$qual_string\n";		
		foreach my $CpG_position(@CpG_positions){
			my $offset = $CpG_position-$read_start+1;
			# print "$target_id\t$CpG_position\t$read_start\t$offset\n";
			next if($offset<0 || $offset >= length($seq));	
                # Depend on different alignmentor(bisreadmapper)
			if($read_strand & 0x10)
            {
				my $qual_score = ord(substr($qual_string,$offset+1,1))-$phred;
				next if($readHapInfo{$read_id}->{$CpG_position}->{"qual"} && $readHapInfo{$read_id}->{$CpG_position}->{"qual"} > $qual_score);
				$readHapInfo{$read_id}->{$CpG_position}->{"base"}=$rcTable{substr($seq,$offset+1,1)};			
				$readHapInfo{$read_id}->{$CpG_position}->{"qual"}=$qual_score;	
			}else{
				my $qual_score = ord(substr($qual_string,$offset,1))-$phred;
				next if($readHapInfo{$read_id}->{$CpG_position}->{"qual"} && $readHapInfo{$read_id}->{$CpG_position}->{"qual"} > $qual_score);
				$readHapInfo{$read_id}->{$CpG_position}->{"base"}=substr($seq,$offset,1);
				$readHapInfo{$read_id}->{$CpG_position}->{"qual"}=$qual_score;

			}	
			#print "$offset\t", $readHapInfo{$read_id}->{$CpG_position}, "\n";
		}	
		$read_start_pos_table{$read_id} = $read_start if(!$read_start_pos_table{$read_id} || $read_start_pos_table{$read_id}<$read_start);
	}		
	
	#Use read start as the UMI to collapse multiple clonal reads.
    # How to define UMI will influence the result compared with classic methylation level(bismark)
	my @read_list = keys(%readHapInfo);		
	my %unique_read_base_info;
	foreach my $read_id (@read_list){
	my $UMI;
	if($library=="RRBS"){
	$UMI = $read_id; # read_id as UMI 
	}else{                        
	$UMI = $read_start_pos_table{$read_id};  # start postion as UMI
	}
	foreach my $CpG_position (keys(%{$readHapInfo{$read_id}})){
		next if(!$readHapInfo{$read_id}->{$CpG_position}->{"base"});
		#print "$CpG_position\n";
		#print $readHapInfo{$read_id}->{$CpG_position}->{"base"}, "\n";
		#print $readHapInfo{$read_id}->{$CpG_position}->{"qual"}, "\n";
		$unique_read_base_info{$UMI}->{$CpG_position}->{$readHapInfo{$read_id}->{$CpG_position}->{"base"}}+=$readHapInfo{$read_id}->{$CpG_position}->{"qual"};
		}
	}
	
	#Derive the consensus haplotype string based on multiple clonal reads.
	foreach my $UMI (keys(%unique_read_base_info)){
		my $hap_string="";
		foreach my $CpG_position(@CpG_positions){
			my $base = "N";
			my $best_qual = 0;
			if($unique_read_base_info{$UMI}->{$CpG_position}){
				foreach my $base_call (keys(%{$unique_read_base_info{$UMI}->{$CpG_position}})){
					next if($unique_read_base_info{$UMI}->{$CpG_position}->{$base_call} <= $best_qual);
					$best_qual = $unique_read_base_info{$UMI}->{$CpG_position}->{$base_call};
					$base = $base_call;
				}
			}					
			$hap_string=$hap_string.$base;
		}	
		$hapInfoTable{$target_id}->{"hap_counts"}->{$hap_string}++ if($hap_string);
	}

	#report haplotype strings only on the positions with valid calls.
	foreach my $hap_string(keys %{$hapInfoTable{$target_id}->{"hap_counts"}}){
		my @valid_positions;
		my $valid_hap="";
		for(my $i=0; $i<length($hap_string); $i++){
			my $allele = substr($hap_string,$i,1);
			next if($allele eq "N");
			$valid_hap = $valid_hap.$allele;
			push(@valid_positions,$CpG_positions[$i]);
		}
		next if(scalar(@valid_positions)<1);
		next if($valid_hap =~ /[AG]/);
		print $target_id,"\t$valid_hap\t", $hapInfoTable{$target_id}->{"hap_counts"}->{$hap_string},"\t",join(",",@valid_positions), "\n";
	}
}
}


   ### Produce Run Time
   my $end_run = time();
   my $run_time = $end_run - $start_run;
   my $days  = int($run_time/(24*60*60));
   my $hours = ($run_time/(60*60))%24;
   my $mins  = ($run_time/60)%60;
   my $secs  = $run_time%60;
   warn "\n# Methtype completed in ${days}d ${hours}h ${mins}m ${secs}s\n";



sub process_command_line{

	my $target_list_file;
	my $merged_bam_file;
	my $library;
	my $Aligner;
	my $chrsizeFile;
	my $cpg_position_file;
	my $genome;
	my $help;
	my $version;
	my $bismark_version;
	
	my $pbs;
	my $submit;
	my %walltime;
	my %ppn;
	my %multicore;
	my $nodes;
	
	my $ppn;
        my $walltime;
        my $multicore;
    
	my $command_line=GetOptions (
		     "bam=s"			=> \$merged_bam_file,
                     "bed=s"   			=> \$target_list_file,
	             "library=s"		=> \$library,
		     "genome=s"			=> \$genome,
		     "cpgpostion=s"		=> \$cpg_position_file,
		     "chromsizes=s"		=> \$chrsizeFile,
		     "aligner=s" 		=> \$Aligner,
	             "phred=s"                  => \$phred,
                     "server"   		=> \$server,
	             "queue"   			=> \$queue,                                                                    
                     "help"      		=> \$help,
                     "pbs"		        => \$pbs,
                     "walltime"			=> \$walltime,
                     "multicore"		=> \$multicore,             
	             "submit"   		=> \$submit,
	             "nodes"   		        => \$nodes, 
		     "version"                  => \$version,
	             );

  #################################################################################################
    ##################### SmartBismark Version and Usage (Version and Usage) ########################
    #################################################################################################
    if ($help){
    print_helpfile();
    exit;
    }  

    if ($version){
    bam2hapinfo_version();
    exit;
    }  

 
    #################################################################################################
    ##################### Assemble Bismark Reference (hg19,hg38,mm9,mm10) ##########################
    #################################################################################################
    my $ip=ipconfig();
    $server="GM" if $ip eq "132.239.25.238";
    $server="TSCC" if $ip eq "132.249.107.88";

    if($server eq "TSCC"){
    	if($genome eq "hg19"){
		$cpg_position_file="/home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt" if !defined $cpg_position_file;
		$chrsizeFile="/home/shg047/oasis/db/hg19/hg19.chrom.sizes" if !defined $chrsizeFile;
		print "Bismark alignment reference: $cpg_position_file\n";
		print "Chromosome Size File: $chrsizeFile\n";
		}elsif($genome eq "hg38"){
		$cpg_position_file="/home/shg047/oasis/db/hg38/HsGenome38.CpG.positions.txt" if !defined $cpg_position_file;
		$chrsizeFile="/home/shg047/oasis/db/hg38/hg38.chrom.sizes" if !defined $chrsizeFile;
		print "Bismark alignment reference: $cpg_position_file\n";
		print "Chromosome Size File: $chrsizeFile\n";
		}elsif($genome eq "mm9"){
		$cpg_position_file="/home/shg047/oasis/db/mm9/mm9.CpG.positions.txt" if !defined $cpg_position_file;
		$chrsizeFile="/home/shg047/oasis/db/mm9/mm9.chrom.sizes" if !defined $chrsizeFile;
		print "Bismark alignment reference: $cpg_position_file\n";
		print "Chromosome Size File: $chrsizeFile\n";
		}elsif($genome eq "mm10"){
		$cpg_position_file="/home/shg047/oasis/db/mm10/mm10.CpG.positions.txt" if !defined $cpg_position_file;
		$chrsizeFile="/home/shg047/oasis/db/mm10/mm10.chrom.sizes" if !defined $chrsizeFile;
		print "Bismark alignment reference: $cpg_position_file\n";
		print "Chromosome Size File: $chrsizeFile\n";
		}else{
		warn("# Please assign genome version (in TSCC)to the script: hg19? hg38? mm9? mm10?");	
		}
    }elsif($server eq "GM"){
		if($genome eq "hg19"){
		$cpg_position_file="/media/Home_Raid1/shg047/work/db/hg19/HsGenome19.CpG.positions.txt" if !defined $cpg_position_file;
		$chrsizeFile="/media/Home_Raid1/shg047/work/db/hg19/hg19.chrom.sizes" if !defined $chrsizeFile;
		print "Bismark alignment reference: $cpg_position_file\n";
		print "Chromosome Size File: $chrsizeFile\n";
		}elsif($genome=="hg38"){
		$cpg_position_file="/media/Home_Raid1/shg047/work/db/hg38/HsGenome38.CpG.positions.txt" if !defined $cpg_position_file;
		$chrsizeFile="/media/Home_Raid1/shg047/work/db/hg38.chrom.sizes" if !defined $chrsizeFile;
		print "Bismark alignment reference: $cpg_position_file\n";
		print "Chromosome Size File: $chrsizeFile\n";
		}elsif($genome eq "mm9"){
		$cpg_position_file="/media/Home_Raid1/shg047/work/db/mm9/mm9.CpG.positions.txt" if !defined $cpg_position_file;
		$chrsizeFile="/media/Home_Raid1/shg047/work/db/mm9/mm9.chrom.sizes" if !defined $chrsizeFile;
		print "Bismark alignment reference: $cpg_position_file\n";
		print "Chromosome Size File: $chrsizeFile\n";
		}elsif($genome eq "mm10"){
		$cpg_position_file="/media/Home_Raid1/shg047/work/db/mm10/mm10.CpG.positions.txt" if !defined $cpg_position_file;
		$chrsizeFile="/media/Home_Raid1/shg047/work/db/mm10/mm10.chrom.sizes" if !defined $chrsizeFile;
		print "Bismark alignment reference: $cpg_position_file\n";
		print "Chromosome Size File: $chrsizeFile\n";
		}else{
		warn("# Please assign genome version (in Genome-miner)to the script: hg19? hg38? mm9? mm10?");	
		}
	}else{
		print "Please provide the cpg_position_file to the script: --cpg_position_file if !defined $cpg_position_file\n";
		print "Please provide the cpg_position_file to the script: --chrsizeFile if !defined $chrsizeFile\n";
	}
    
    #################################################################################################
    ##################### Assemble PBS Paramters (nodes, ppn, walltime multicore) ##########################
    #################################################################################################

    unless (defined $merged_bam_file and defined $target_list_file and defined $Aligner and defined $library and defined $chrsizeFile and defined $library){
    warn "\n\n\t# Error: Please respecify essential command line options!\n";
    print_helpfile();
    }
    #################################################################################################
    ##################### Assemble PBS Paramters (nodes, ppn, walltime multicore) ##########################
    #################################################################################################
    warn "\n# You didn't assign --queue for the script, the default setting: multicore=2 and ppn=6 will be applied!\n\n" if ! defined $queue; 

    $queue="hotel" if ! defined $queue;
    %walltime=(
    hotel   => "168:00:00",
    condo   => "8:00:00",
    pdafm   => "72:00:00",
    glean   => "72:00:00",
    default => "168:00:00",
    );
    %ppn=(
    hotel   => "16",
    pdafm   => "32",
    glean   => "16",
    condo   => "16",
    default => "6",
    );
    %multicore=(
    hotel   => "6",
    pdafm   => "12",
    glean   => "6",
    condo   => "6",
    default => "2",
    );
    warn "# Queue: $queue is not found in this server, please check the queue name\n" if ! defined $ppn{$queue}; 
    $nodes=1;
    $ppn=$ppn{$queue};
    $walltime=$walltime{$queue};
    $multicore=$multicore{$queue};
  
    #################################################################################################
    ##################### Return All the paramters for SmartBismark ##########################
    #################################################################################################
    return($merged_bam_file,$target_list_file,$library,$cpg_position_file,
           $chrsizeFile,$Aligner,$phred,$server,$queue,$help,$pbs,$submit,$nodes,$ppn,
           $walltime,$multicore);
           
}


sub bam2hapinfo_version{
	
	print<<"VERSION";
	
          bam2hapinfo - Smart tools to extract methylation haplotypes from BAM files (BisReadMapper and Bismark).

                               Bam2hapinfo Version: "bam2hapinfo version: 1.02"

			            BisReadMapper : <hdinhdp\@gmail.com>
                           
                 Copyright 2010-15 Kun Zhang <kzhang\@eng.ucsd.edu >, University of California, San Diego

			           Software Maintainer: <shg047\@ucsd.edu>

VERSION
exit;
}

sub pbsprocess($$$$$$$$){
	my $merged_bam_file=${shift @_};
	my $target_list_file=${shift @_};
	my $library=${shift @_};
	my $cpg_position_file=${shift @_};
	my $chrsizeFile=${shift @_};
	my $Aligner=${shift @_};
	my $phred=${shift @_};
	my $queue=${shift @_};
	my $ppn=${shift @_};
	my $nodes=${shift @_};
	my $walltime=${shift @_};	
	my $dir=${shift @_};
	
	print "#!/bin/csh\n";
    print "#PBS -N $merged_bam_file\n";
    print "#PBS -q $queue\n";  # glean is free
    print "#PBS -l nodes=$nodes:ppn=$ppn\n";
    print "#PBS -l walltime=$walltime\n";
    print "#PBS -o $merged_bam_file.log\n";
    print "#PBS -e $merged_bam_file.err\n";
    print "#PBS -V\n";
    print "#PBS -m abe\n";
    print "cd $dir\n"; 
    print "perl ~/bin/bam2hapinfo.pl --bam $merged_bam_file --bed $target_list_file --library $library --aligner $Aligner --phred $phred --cpgfile $cpg_position_file --chrsizeFile $chrsizeFile\n";
	exit;
}


sub bismark_version{
	my @version=`bismark --version`;
	my $bismark_version;
	foreach my $line(@version){
		if($line=~/version/i){
			$bismark_version=$line;
		}
	}
	return($bismark_version);
}

sub ipconfig{
my $sock = IO::Socket::INET->new(
                       PeerAddr=> "example.com",
                       PeerPort=> 80,
                       Proto   => "tcp");
my $localip = $sock->sockhost;
return($localip)
}

sub print_helpfile{

 print<<"HOW_TO";

DESCRIPTION

USAGE: bam2hapinfo --bam Indx01.bam --bed mhb.bed --library RRBS --aligner bismark --phred 33 --cpgfile HsGenome19.CpG.positions.txt --chrsizeFile hg19.chrom.sizes 

Last edited on 15 Apirl 2017 Shicheng Guo <Shihcheng.Guo\@gmail.com>.
	
ARGUMENTS:
		
 --bam              Input bam files aligned by Bismark, BisReadMapper. We will updated more aligners later.
	
 --bed              Genomic regions of interesting, such as MHB, Enhancer or your custom defined regions.

 --library          Defined methylation sequencing library types: WGBS(no capture based), RRBS(capture based)

 --aligner          Assign alignment tools of the BAM files, Now support: Bismark, BisReadMapper.

 --phred            Assign phred score for BAM files. Phred=33 or 64. Default setting is phred=33;
 
 --genome           Reference genome: mm9,mm10,hg19 and hg38

 --cpgfile          CpG positions for each genome (mm9,mm10,hg19 and hg38).

 --chrsizeFile      chrsizeFile can be download with UCSC tools like:  fetchChromSizes hg19 > hg19.chrom.sizes


Other options:   

 --server           TSCC, GM(GenomeMiner). Combined with genome, --server and --genome will provide the 
	            location of the alignment reference for bismark.  
	
 --submit           Submit pbs job or not. SmartBismark will creat pbs job files for each fastq file and defaulty
                    taken the system is PBS system. if --submit="submit", then PBS job will be submitted and PBS 
                    ID will be printed in the STANDOUT.    

 --queue            TSCC queue: hotel, glean, pdafm, condo		
 
 --nodes            TSCC nodes: nodes=3,6,9 and maximum 24.		
 
    
	
HOW_TO
exit;
}


	

