

#!/usr/bin/perl
use strict;
use Cwd;

chdir getcwd;

my $file="WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.time3.upstream.bed ";
my @ref=qw(LOCK.hg19.bed 
Coloncancer.small.DMR.PMC3145050.hg19.bed 
VMR.charm.PMID20844285.hg19.bed 
Hic.topological.domain.hESC.hg19.bed 
Hic.boundary.IMR90.hg19.bed 
Hic.common.boundary.hESC.IMR90.hg19.bed 
Hic.topological.domain.IMR90.hg19.bed
)


foreach(my $i in @ref){
my $a=`bedtools intersect -wa -u  -a $file -b ../../annotation/hg19/$i | wc -l`;
my $b=`wc -l $file`;

my $c=`bedtools intersect -wa -u  -b $file -a ../../annotation/hg19/$i | wc -l`;
my $d=`wc -l $i`;

my $e1=`bedtools intersect -wa -u  -b $file -a ../../annotation/hg19/$i | wc -l`;

chomp($a);
print "$i\t$a\n";
}




bedtools intersect -wa -u  -a  WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.time3.upstream.bed -b ../../annotation/hg19/LOCK.hg19.bed | wc -l
bedtools intersect -wa -u  -a  WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.time3.upstream.bed -b ../../annotation/hg19/Coloncancer.small.DMR.PMC3145050.hg19.bed | wc -l
bedtools intersect -wa -u  -a  WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.time3.upstream.bed -b ../../annotation/hg19/VMR.charm.PMID20844285.hg19.bed | wc -l
bedtools intersect -wa -u  -a  WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.time3.upstream.bed -b ../../annotation/hg19/Hic.topological.domain.hESC.hg19.bed | wc -l
bedtools intersect -wa -u  -a  WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.time3.upstream.bed -b ../../annotation/hg19/Hic.boundary.IMR90.hg19.bed | wc -l
bedtools intersect -wa -u  -a  WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.time3.upstream.bed -b ../../annotation/hg19/Hic.boundary.hESC.hg19.bed | wc -l
bedtools intersect -wa -u  -a  WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.time3.upstream.bed -b ../../annotation/hg19/Hic.common.boundary.hESC.IMR90.hg19.bed | wc -l
bedtools intersect -wa -u  -a  WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.time3.upstream.bed -b ../../annotation/hg19/Hic.topological.domain.IMR90.hg19.bed | wc -l
