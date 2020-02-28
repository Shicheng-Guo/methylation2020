#!/usr/bin/perl
use strict;
use Cwd;
die &USAGE if scalar(@ARGV)<1;
my $bam_dir=getcwd;
my $bed_dir="/home/shg047/oasis/mouse/RD/Mouse.MHB.Alice.chrY.bam.RD90_80up.bed";
my $filter=shift @ARGV;
chdir $bam_dir;
my @file=glob("*.bam");
foreach my $file(@file){
my ($sample,undef)=split /$filter/,$file;       # for monod dataset
my $chr=$1 if $file=~/(chr\w+)/;
print "$chr\n";
#my ($sample,undef)=split /_/,$file;            # for Dennis Lo's dataset
#my ($sample,undef)=split /.sorted.clipped.bam/,$file;            # for Dennis Lo's dataset
print "$sample\t$bam_dir/$file\t$bed_dir\n";
}

sub USAGE{
print "USAGE: perl $0 [Bam_Directory]\n";
print '
Sample_Info_File Format: 
ENCFF000LUN_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LUN_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LUP_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LUP_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LUQ_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LUQ_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LUT_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LUT_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LUU_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LUU_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LUV_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LUV_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LVA_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LVA_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LVB_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LVB_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LVE_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LVE_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LVF_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LVF_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LVI_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LVI_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LVJ_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LVJ_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
'
}
