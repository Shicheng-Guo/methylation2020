#!/usr/bin/perl -w

# MasterBisReadMapper.pl: a perl script to process bisulfite reads or mapped sam files
# Contact: Dinh Diep
# Version 1.4

use strict;
use warnings;
use Getopt::Std;

my %read_files;
my %sam_files;
my %bam_files;
my %chroms;

my ($scripts_dir, $samtools, $bamutils, $samtools_snp) = (0,0,0,0);
my ($cpg_list, $snp_list, $target_bed) = (0,0,0);
my ($ref_dbsnp, $ref_fai, $ref_fa) = (0,0,0);
my ($mapper, $trimgalore) = (0,0);
my ($template_fwd, $template_rev) = (0,0);

sub main{
	my %opts = ();
	getopts('i:s:v:b:p:d:c:m:', \%opts);
	die( printUsage() ) if(@ARGV == 0 and !$opts{i});
	die("[MasterBisReadMapper] No sam file provided.\n") if(!$opts{i});
	open(INFILE, "$opts{i}") || die("Cannot open input file(s) list\n");
	while(my $line = <INFILE>){
		chomp($line);
		my @tmp = split /\t/, $line;
		my ($dir, $file) = ($tmp[1], $tmp[2]);
		if($file =~ m/sam$/){
			# file_path sample_id base_quality
			if($opts{c} and $opts{c} eq "yes"){
				open(CHR_SAM, "$dir/$file") || die("[MasterBisReadMapper] Error reading $dir/$file\n");
				my $tmp_line = <CHR_SAM>;
				close(CHR_SAM);
				my @fields = split /\t/, $tmp_line;
				my $cur_chr = $fields[2];
				push(@{$sam_files{$tmp[0].".".$cur_chr}->{"data"}}, $dir."\t".$file."\t".$tmp[3]."\t".$tmp[4]);
				
			}else{
				push(@{$sam_files{$tmp[0]}->{"data"}}, $dir."\t".$file."\t".$tmp[3]."\t".$tmp[4]);
			}
		}else{
			if(!$tmp[7]){ die("No trim mode specified for $file !!!!\n"); }
			# <sample id> <dir> <read1.fq | read1.fq,read2.fq | *.sam> <phred> <clonal method> <adaptor r1> <adaptor r2> <trim mode>
			push(@{$read_files{$tmp[0]}->{"data"}}, $dir."\t".$file."\t".$tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[6]."\t".$tmp[7]);
		}
	}
	close(INFILE);
	die("[MasterBisReadMapper] No path file provided.\n") if(!$opts{s});
	open(INFILE, "$opts{s}") || die("[MasterBisReadMapper] Cannot open path file \n");
	while(my $line = <INFILE>){
		chomp($line);
		my @f = split /[="]/, $line;
		next if(!$f[0]);
		$scripts_dir=$f[2] if($f[0] eq "scripts_dir");
		$samtools=$f[2] if($f[0] eq "samtools");
		$bamutils=$f[2] if($f[0] eq "bamUtils");
		$samtools_snp=$f[2] if($f[0] eq "samtools_snp");
		$cpg_list=$f[2] if($f[0] eq "cpg_list");
		$snp_list=$f[2] if($f[0] eq "snp_list");
		$ref_dbsnp=$f[2] if($f[0] eq "ref_dbsnp");
		$ref_fai=$f[2] if($f[0] eq "ref_fai");
		$ref_fa=$f[2] if($f[0] eq "ref_fa");
		$template_fwd=$f[2] if($f[0] eq "template_fwd");
		$template_rev=$f[2] if($f[0] eq "template_rev");
		$target_bed=$f[2] if($f[0] eq "target_bed");
		$mapper=$f[2] if($f[0] eq "mapper");
		$trimgalore=$f[2] if($f[0] eq "trimgalore");
	}
	close(INFILE);
	die("[MasterBisReadMapper] No samtools_snp path but -v is specified.\n") if($opts{v} and $opts{v} eq "yes" and !$samtools_snp);
	die("[MasterBisReadMapper] No samtools path.\n") if(!$samtools);
	die("[MasterBisReadMapper] ref_fai, or ref_fa missing.\n") if(!$ref_fai || !$ref_fa );
	die("[MasterBisReadMapper] No CpG position list given.\n") if(!$cpg_list);
	
	#### Process fastq files #### 
	foreach my $name (keys %read_files){
		my $out_dir = $name . "_map";
		my $cmd = "mkdir $out_dir";
		system($cmd) == 0 or warn "[MasterBisReadMapper] Failed to make a new directory to store files (exit $?): $!\nCommand used:\n\t$cmd\n";
		my @records = @{$read_files{$name}->{"data"}};
		for(my $i = 0; $i < scalar(@records); $i++){
			my ($dir, $read1_and_read2, $phred_base, $clonal_method, $adaptor_r1, $adaptor_r2, $trim_mode) = split /\t/, $records[$i];
			$trim_mode = uc($trim_mode);
			my ($read1, $read2) = ($read1_and_read2, 0);
			($read1, $read2) = split ",", $read1_and_read2 if($read1_and_read2 =~ s/,/,/g);
			my $fqname = $read1;
			$fqname =~ s/.fq//g;
			$fqname =~ s/.fastq//g;
			$fqname = $fqname. ".PE" if($read2);

			my $err_status = process_fastq($out_dir, $dir, $fqname, $read1, $read2, $phred_base, $clonal_method, $adaptor_r1, $adaptor_r2, $trim_mode);
		        my @list_sams;
			if($err_status){
				print "[MasterBisReadMapper] !!!! Error while processing $read1_and_read2 in $dir !!!!\n";
				next;
			}
			my $clonal_id_method = $clonal_method;
			$clonal_id_method = "umi" if($clonal_method =~ m/umi/i);
			print "[MasterBisReadMapper] Here are the sam files generated:\n";
			open(IN_STATUS, "$out_dir/$fqname.status") || warn("[MasterBisReadMapper] Cannot read status file: $out_dir/$fqname.status\n");
			while(my $line = <IN_STATUS>){
				chomp($line);
				if($line =~ m/^chr/){
					my @tmp = split "\t", $line;
					my $chr = $tmp[0];
					my $file = "$fqname.$chr.sorted.sam";
					if($opts{c} and $opts{c} eq "yes"){
						my $cur_name = $name ."." . $chr;
						push(@{$sam_files{$cur_name}->{"data"}}, $out_dir."\t".$file."\t33\t".$clonal_id_method);
						print "$cur_name\t$out_dir\t$file\t33\t$clonal_method\n";
					}else{
						push(@{$sam_files{$name}->{"data"}}, $out_dir."\t".$file."\t33\t".$clonal_id_method);
						print "$name\t$out_dir\t$file\t33\t$clonal_method\n";
					}
				}
			}
			close(IN_STATUS);
		}
	}
	print "[MasterBisReadMapper] All reads have been mapped \n";
	exit 0 if($opts{m} and $opts{m} eq "yes");
	#### Process sam files ####
	foreach my $name (keys %sam_files){
		my $cmd = "mkdir -p $name"."_map/methylfiles";
		system($cmd) == 0 or warn "[MasterBisReadMapper] Failed to make a new directory to store files (exit $?): $!\nCommand used:\n\t$cmd\n";
		my $depth = $opts{d};
		$depth = 5 if(!$depth);
		my $pileup_opt = $opts{p};
		$pileup_opt = "no" if(!$pileup_opt);
		my $variant_opt = $opts{v};
		$variant_opt = "no" if(!$variant_opt);
		my $bedfile_opt = $opts{b};
		$bedfile_opt = "no" if(!$bedfile_opt);
		my $err_status = sam_to_methyl("$name"."_map/methylfiles/", $name, $depth, $pileup_opt, $variant_opt, $bedfile_opt);
		if($err_status){
			print "[MasterBisReadMapper] !!!! Error while processing SAM files for $name !!!! \n";
		}
	}
	
}

sub process_fastq{
	my ($out_dir, $dir, $name, $read1, $read2, $phred_base, $clonal_method, $adaptor_r1, $adaptor_r2, $trim_mode) = @_;
	my @list_sams;
	my $cmd;
	my $trimOpts;
	my $trimOpts_wgbs_se = "--clip_R1 5";
	my $trimOpts_rrbs_se = "--rrbs --non_directional";
	my $trimOpts_bspp_se = "--clip_R1 27";
        my $trimOpts_wgbs_pe = "--paired --clip_R1 5 --clip_R2 5";
        my $trimOpts_rrbs_pe = "--paired --rrbs --non_directional";
        my $trimOpts_bspp_pe = "--paired --clip_R1 27 --clip_R2 27";
	if($read2){
		$trimOpts = $trimOpts_wgbs_pe if($trim_mode eq "WGBS");
		$trimOpts = $trimOpts_rrbs_pe if($trim_mode eq "RRBS");
		$trimOpts = $trimOpts_bspp_pe if($trim_mode eq "BSPP");
	}else{	
		$trimOpts = $trimOpts_wgbs_se if($trim_mode eq "WGBS");
		$trimOpts = $trimOpts_rrbs_se if($trim_mode eq "RRBS");
		$trimOpts = $trimOpts_bspp_se if($trim_mode eq "BSPP");
	}
	#### get UMI ####
	if($clonal_method =~ m/UMI/i){
		my ($tmp, $length) = split /\:/, $clonal_method;
		print "[MasterBisReadMapper] $name, UMI length : $length\n";
		if($read2){
			$cmd = "$scripts_dir/getUMI.pl $dir/$read1 $dir/$read2 $out_dir/$name $length";
			if(system($cmd) != 0){warn "[MasterBisReadMapper] Failed to get UMIs (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
			$read1 = $name.".R1.fq";
			$read2 = $name.".R2.fq";
			$dir = $out_dir;
		}else{
			$cmd = "$scripts_dir/getUMI.pl $dir/$read1 $out_dir/$name $length";
			if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to get UMIs (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
			$read1 = $name.".R1.fq";
			$dir = $out_dir;
		}
	}
	
	#### trim galore ####
	if($read2){
		$cmd = "$trimgalore --phred$phred_base -o $out_dir --suppress_warn --dont_gzip $trimOpts -a $adaptor_r1 -a2 $adaptor_r2 $dir/$read1 $dir/$read2";
		if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to run trimgalore (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
		#unlink("$dir/$read1");
		#unlink("$dir/$read2");
		$read1 =~ s/\.fq$//g;
		$read1 =~ s/\.fastq$//g;
		$read1 =~ s/\.fq.gz$//g;
		$read1 =~ s/\.fastq.gz$//g;
		
		$read1 = $read1."_val_1.fq";
		$read2 =~ s/\.fq$//g;
		$read2 =~ s/\.fastq$//g;
		$read2 =~ s/\.fq.gz$//g;
		$read2 =~ s/\.fastq.gz$//g;
		$read2 = $read2."_val_2.fq";
		$dir = $out_dir;
		$cmd = "cat $dir/$read2 >> $dir/$read1";
		if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to concatenate all read files (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
		#unlink("$dir/$read2");
	}else{
		$cmd = "$trimgalore --phred$phred_base -o $out_dir --suppress_warn --dont_gzip $trimOpts -a $adaptor_r1 $dir/$read1";
		if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to run trimgalore (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
		#unlink("$dir/$read1");
		$read1 =~ s/\.fq$//g;
                $read1 =~ s/\.fastq$//g;
		$read1 = $read1."_trimmed.fq";
		$dir = $out_dir;
	}
	
	#### bisreadmapper ####
	$cmd = "$scripts_dir/BisReadMapper.pl -r $dir/$read1 -W $template_fwd -C $template_rev -g $ref_fai -a $mapper -b $phred_base -p 4 -n $out_dir/$name > $out_dir/$name.status";
	print "[MasterBisReadMapper] Begins BisReadMapper...\n";
	if(system($cmd) != 0){warn "[MasterBisReadMapper] Failed to run BisReadMapper (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
	#unlink("$dir/$read1");
	return 0;
}

sub sam_to_methyl{
	my $out_dir = shift;
	my $name = shift;
	my $depth = shift;
	my $pileup_opt = shift;
	my $variant_opt = shift;
	my $bedfile_opt = shift;
	my $bam_file = $out_dir . "/" . $name;
	my $sam_file = $out_dir . "/" . $name . ".merged.sam";
	unlink($sam_file);
	my $pileup_file = $out_dir . "/" . $name . ".pileup";
	my $methylFreq_file = $out_dir . "/" . $name . ".methylFreq";
	my $bed_file = $out_dir ."/". $name . ".BED.txt";
	my $snp_file = $out_dir ."/". $name . ".SNP.txt";
	my $filtered_snp_file = $out_dir. "/". $name . ".filtered.SNP.txt";
	my ($dir, $file, $phred_base, $clonal_id_method) = ("NA", "NA", 33, "none");
	# ==== combine all the sam files to one ===== #
	foreach my $record (@{$sam_files{$name}->{"data"}}){
		chomp($record);
		print "$name\t$record\n";
		($dir, $file, $phred_base, $clonal_id_method) = split /\t/, $record;
		#process this file
		my $cmd = "cat $dir/$file >> $sam_file";
		if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to merge sam files (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
	}	
	# ==== convert to bam, sort, and rmdup ==== #
	my $cmd = "$samtools view -ubSt $ref_fai $sam_file > $bam_file.bam";
	if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to convert bam file (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
	unlink("$sam_file");
	
	# ==== fix mate information ===== #
	$cmd = "$samtools view $bam_file.bam | sort -k 1,1 | $scripts_dir/fixMateInfo.pl | $samtools view -bSt $ref_fai - | $samtools sort - $bam_file.sorted";
	#$cmd = "$samtools sort $bam_file.bam $bam_file.sorted";
	if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to fixMateInfo (exit $?): $!\nCommand used:\n\t$cmd\n";return 1;}
	unlink("$bam_file.bam");

	$bam_file = $bam_file.".sorted";

	print "[MasterBisReadMapper] $name mapped reads stats\n";
	$cmd = "$samtools flagstat $bam_file.bam";
	if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to run flagstats (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}

	# ==== clip overlaps ==== #
	print "[MasterBisReadMapper] $name clipping PE reads overlaps\n";
	$cmd = "$bamutils clipOverlap --stats --poolSize 4000000 --in $bam_file.bam --out $bam_file.clipped.bam";
	#$cmd = "mv $bam_file.bam $bam_file.clipped.bam";
	if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to clipOverlap (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
	unlink("$bam_file.bam");
	$bam_file = $bam_file . ".clipped";

	# ==== remove duplicates ===== #
	if($clonal_id_method eq "umi"){
		print "[MasterBisReadMapper] $name after duplicates removed:\n";
		$cmd = "$samtools view $bam_file.bam | $scripts_dir/sam_UMI_filter_Poisson.pl | $samtools view -bSt $ref_fai - > $bam_file.rmdup.bam";
		if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to remove duplicates with UMI filter (exit $?): $!\nCommand used:\n\t$cmd\n";return 1;}
		unlink("$bam_file.bam");
		$bam_file = $bam_file . ".rmdup";
		$cmd = "$samtools flagstat $bam_file.bam";
		if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to run flagstats (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
	}
	if($clonal_id_method eq "samtools"){
		print "[MasterBisReadMapper] $name after duplicates removed:\n";
		$cmd = "$samtools rmdup -S $bam_file.bam $bam_file.rmdup.bam";
		if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to remove duplicates with samtools (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
		unlink("$bam_file.bam");
		$bam_file = $bam_file . ".rmdup";
		
		$cmd = "$samtools flagstat $bam_file.bam";
		if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to run flagstats (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
	}
	$bam_file = $bam_file .".bam";
	$cmd = "$samtools index $bam_file";
	if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to index bam file (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
	# ==== extract pileup, methylFreq, bed, and snps === #
	if($pileup_opt eq "yes"){
		$cmd = "$samtools mpileup -BA -f $ref_fa $bam_file > $pileup_file";
		if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to generate pileup (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
		$cmd = "$scripts_dir/extractMethyl.pl $cpg_list $phred_base < $pileup_file > $methylFreq_file";
		if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to extract methyl from pileup (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}

		# calculate on target here
		if($target_bed and $target_bed ne "NA"){
			$cmd = "$scripts_dir/getProbeBias.pl $target_bed < $pileup_file > $methylFreq_file";
			if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to calculate on target bases (exit $?): $!\nCommand used:\n$cmd\n";}
		}

	}else{
		$cmd = "$samtools mpileup -BA -f $ref_fa $bam_file | $scripts_dir/extractMethyl.pl $cpg_list $phred_base > $methylFreq_file";
		if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to extract methyl from bam (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
		
		# calculate on target here
		if($target_bed and $target_bed ne "NA"){
                        $cmd = "$samtools mpileup -BA -f $ref_fa $bam_file | $scripts_dir/getProbeBias.pl $target_bed > $methylFreq_file";
                        if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to calculated on target bases (exit $?): $!\nCommand used:\n\t$cmd\n";}
                }

	}
	$cmd = "$scripts_dir/frMethylCorr.pl $depth < $methylFreq_file";
	system($cmd) == 0 or warn "[MasterBisReadMapper] Failed to perform correlation (exit $?): $!\nCommand used\n\t$cmd\n";

	if($bedfile_opt eq "yes"){
		$cmd = "$scripts_dir/methylFreq2BED.pl $name $depth < $methylFreq_file > $bed_file";
		system($cmd) == 0 or warn "[MasterBisReadMapper] Failed to generate BED (exit $?): $!\nCommand used:\n\t$cmd\n";
	}

	return 0 if($variant_opt eq "no");
	# make separate BAM files
	$cmd = "$samtools view $bam_file | awk '{if(\$2 == 0) print \$0;}' - | $samtools view -uSbt $ref_fai - > $bam_file.Watson";
	if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to generate Watson bam file (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
	$cmd = "$samtools view $bam_file | awk '{if(\$2 == 16) print \$0;}' - | $samtools view -uSbt $ref_fai - > $bam_file.Crick";
	if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to generate Crick bam file (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
	# call the variants
	$cmd = "$samtools_snp pileup -Ac -f $ref_fa $bam_file.Watson | $scripts_dir/extractSNPs.pl W VAR $phred_base > $snp_file";
	if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to call variants on Watson (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
	$cmd = "$samtools_snp pileup -Ac -f $ref_fa $bam_file.Crick | $scripts_dir/extractSNPs.pl C VAR $phred_base >> $snp_file";
	if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to call variants on Crick (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}

	if($snp_list and $snp_list ne "NA"){
		# call homozygous reference (Need a snp list!)
		$cmd = "$samtools_snp pileup -Ac -l $snp_list -f $ref_fa $bam_file.Watson | $scripts_dir/extractSNPs.pl W REF $phred_base >> $snp_file";
		if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to call homozygous reference on Watson (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
		$cmd = "$samtools_snp pileup -Ac -l $snp_list -f $ref_fa $bam_file.Crick | $scripts_dir/extractSNPs.pl C REF $phred_base >> $snp_file";
		if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to call homozygous reference on Crick (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
	}

	unlink("$bam_file.Watson");
	unlink("$bam_file.Crick");
	return 0 if(!$ref_dbsnp);
	# filter snps
	$cmd = "$scripts_dir/bisSnpFilter_DD.pl $scripts_dir $snp_file $ref_dbsnp > $filtered_snp_file";
	if(system($cmd) != 0) {warn "[MasterBisReadMapper] Failed to filter genotypes (exit $?): $!\nCommand used:\n\t$cmd\n"; return 1;}
}

sub printUsage{
	print "MasterBisReadMapper.pl: a perl script to process bisulfite mapped sam files\n";
	print "        -i <list_files>    : <list_sam_files> is a table (tab separated values) of all the fastq and sam files to be processes[Required]\n";
	print "                                 <sample id> <dir> <read1.fq | read1.fq,read2.fq | *.sam> <phred> <clonal method> <adaptor r1> <adaptor r2> <trim mode>\n";
	print "                                 Note: trim mode can be: RRBS, BSPP, or WGBS\n";
	print "                                 Note: umi in sam files are indicated in the reads name field, separated by '_' and is only part with no [0-9]\n";
	print "                                 Note: umi are first X number of bases in read 1, use UMI:X for clonal method, where X is length of UMI.\n";
        print "        -s <list_paths>        : <list_paths> is a file that lists all of the required paths[Required]\n";
	print "        -c [yes/no]            : split processed files by chromosomes or not[Default no]\n";
	print "        -v [yes/no]            : indicate whether to call SNPs or not[Default no]\n";
	print "        -b [yes/no]            : indicate whether to generate BED format file for methylation frequencies[Default no]\n";
	print "        -p [yes/no]            : keep pileup yes or no[Default no]\n";
	print "        -m [yes/no]            : Map only, do not make methylFreq [Default no]\n";
	print "        -d <mindepth>          : minimum depth for BED and frCorr[Default 5]\n";
}

main();
