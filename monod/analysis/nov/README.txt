			    DMAP test data
			    --------------

By now you will have unpacked these files using something like:

gzip -dc mds_test_data.tar.gz | tar -xvf -

producing the directory (this) mds_test_data.  The fastq files are a
restricted set of real patient and control data selected for mapping
to the first 10Mb of human chromosome 21.  The restriction is to
enable fast mapping and checking of the SW and to illustrate the use
of the various DMAP programs.  The files ctrl_*.fastq are from
'normal' individuals and the set mds_*.fastq are from patients with
myelodisplastic syndrome.  Note that these are given only as examples
of program operation and do not imply that any significant MDS-related
methylation changes exist in this region of chromosome 21.

For simplicity, I shall assume that fasta files for the reference
human genome and the bismark/bowtie indices are at
/HomoSapiens/hs_ref_GRCh37/ and that the result files are to be put in
this working directory.  It is also assumed that executables for the
DMAP programs, bismark and bowtie are in the PATH and will execute by
name.  This would be the expected behaviour if these have all been put
in /usr/local/bin in the usual Unix/Linux manner.

1. Mapping: the following commands will map each of the fastq files
against the reference bisulphite genome:

bismark -n 1 /HomoSapiens/hs_ref_GRCh37/ ctrl_1.fastq
bismark -n 1 /HomoSapiens/hs_ref_GRCh37/ ctrl_2.fastq
bismark -n 1 /HomoSapiens/hs_ref_GRCh37/ ctrl_3.fastq
bismark -n 1 /HomoSapiens/hs_ref_GRCh37/ mds_1.fastq
bismark -n 1 /HomoSapiens/hs_ref_GRCh37/ mds_2.fastq
bismark -n 1 /HomoSapiens/hs_ref_GRCh37/ mds_3.fastq

producing a series of files:

ctrl_1.fastq_bismark.sam
ctrl_1.fastq_bismark_SE.alignment_overview.png
ctrl_1.fastq_bismark_SE_report.txt

for each.  On a reasonably well resourced desktop computer, this
should take about 3-4 minutes for each of these fastq files.

This directory contains example sam files from my mapping with names
tst_ctrl_1.fasta_bismark.sam, etc.  It is possible that bismark will
return data in a different order hence simple comparisons of other
runs with these files may show differences.

2. Identifying differentially methylated regions (DMRs), pairwise: the
following command applies Fisher's Exact statistic (-P 40,220) to a
pair of sam files from the above.  The options given require that at
least 2 CpGs in each fragment have 10 or more hits (-F 2 -t 10) and
that the leading CpG of 3' mapped reads is assigned to the previous
fragment (-N).  Only chromosome 21 is to be used (-c 21):

diffmeth -g /HomoSapiens_genome/hs_ref_GRCh37/Homo_sapiens.GRCh37.65.dna.chromosome. \
 -c 21 -P 40,220 -F 2 -t 10 -N -R ctrl_1.fastq_bismark.sam -R mds_1.fastq_bismark.sam

3. Identifying DMRs - ChiSQ test: on a cohort of 6 individuals.  The
following command runs diffmeth using the ChiSQ statistic (-X 40,220)
requiring that at least 2 CpGs in each fragment have 10 or more hits
(-F 2 -t 10), that at least 4 individuals contribute to the statistic
(-I 4) and that the leading CpG of 3' mapped reads is assigned to the
previous fragment (-N).  Only chromosome 21 is to be used (-c 21).
The comparison including control and MDS individuals is intended to
illustrate the use of this statistic with the test data set, not to
imply that there is no difference between the groups:

diffmeth -g /HomoSapiens_genome/hs_ref_GRCh37/Homo_sapiens.GRCh37.65.dna.chromosome. \
 -c 21 -X 40,220 -F 2 -t 10 -N -I 4 -R ctrl_1.fastq_bismark.sam \
 -R ctrl_2.fastq_bismark.sam -R ctrl_3.fastq_bismark.sam \
 -R mds_1.fastq_bismark.sam -R mds_2.fastq_bismark.sam \
 -R mds_3.fastq_bismark.sam

4. Comparing methylation between two groups: this applies the ANOVA F
ratio test (-a 40,220), requiring that at least 4 individuals show
counts for a fragment to be included (-I 4) and the leading CpG of 3'
mapped reads is assigned to the preceding fragment (-N).  Sam data for
the treatment/disease group is identified with -S:

diffmeth -g /HomoSapiens_genome/hs_ref_GRCh37/Homo_sapiens.GRCh37.65.dna.chromosome. \
-c 21 -a 40,220 -N -I 4 -R ctrl_1.fastq_bismark.sam \
-R ctrl_2.fastq_bismark.sam -R ctrl_3.fastq_bismark.sam \
-S mds_1.fastq_bismark.sam -S mds_2.fastq_bismark.sam \
-S mds_3.fastq_bismark.sam

5. As for 4, but indicate the more methylated group (-A 40,220) and
restrict the output to Pr < 0.01 (-U 0.01).  Each line is suffixed
with a 'R' or 'S' character to indicate which group had higher
methylation and a summary of the valid sample counts for R & S groups.

diffmeth -g /HomoSapiens_genome/hs_ref_GRCh37/Homo_sapiens.GRCh37.65.dna.chromosome. \
-c 21 -A 40,220 -U 0.01 -N -I 4 -R ctrl_1.fastq_bismark.sam \
-R ctrl_2.fastq_bismark.sam -R ctrl_3.fastq_bismark.sam \
-S mds_1.fastq_bismark.sam -S mds_2.fastq_bismark.sam \
-S mds_3.fastq_bismark.sam

----------------------------------------------------------------------
Peter Stockwell
27-Jan=2014
