The following are usage instructions for the three main programs:

Step-by-step:

1) use genomePrep.pl to generate the in-silico bisulfite converted references, C->T for Watson, and G->A for Crick strands. Note that bisulfite conversion makes the two strands non-complementary.

2) use aligner software to create index files from the reference sequence for alignment. Supported aligners are bowtie2 (version bowtie2-2.1.0), bwa (bwa-0.7.5a), SOAP2 (soap2.21release) , LAST (last-458) , or GEM (GEM-binaries-Linux-x86_64-core_i3-20130406-045632)
** NOTE: For reads longer than 100 bp, we recommend using bwa mem
** NOTE: Both strands (*.bis.CT and *.bis.GA) can be concatenated into one file and only one index needs to be created so long as the aligner can support larger index files.

3) use samtools to generate the *.fai file from the reference sequence.

4) generate a table of the files to be processed
<sample id> <dir> <read1.fq | read1.fq,read2.fq | *.sam> <phred> <clonal method> <adaptor r1> <adaptor r2>	

5) generate the list_paths file

6) run MasterBisReadMapper.pl with list_sam_files and list_paths as inputs.


===================================================================
Usage: genomePrep.pl genome.fa[.gz] context=[cg/all] convert=[yes/no]
convert: yes = convert genome, no = don't convert genome
context: CG = CG only, ALL = all C context
** This program will generate the CpG position files when the context option is set to CG or ALL
===================================================================
Usage: MasterBisReadMapper.pl -i </path/to/list_sam_file> -s </path/to/list_path_file> [options] &> log
Required:
	-i <list_files>    : <list_sam_files> is a table (tab separated values) of all the fastq and sam files to be processes[Required]
                                 <sample id> <dir> <read1.fq | read1.fq,read2.fq | *.sam> <phred> <clonal method> <adaptor r1> <adaptor r2> <TRIM mode>
                                 Note: umi in sam files are indicated in the reads name field, separated by '_' and is only part with no [0-9]
                                 Note: umi are first X number of bases in read 1, use UMI:X for clonal method, where X is length of UMI.
        -s <list_paths>        : <list_paths> is a file that lists all of the required paths[Required]
Optional settings:
	-c [yes/no]            : split processed files by chromosomes or not[Default no]
        -v [yes/no]            : indicate whether to call SNPs or not[Default no]
        -b [yes/no]            : indicate whether to generate BED format file for methylation frequencies[Default no]
        -p [yes/no]            : keep pileup yes or no[Default no]
        -m [yes/no]            : Map only, do not make methylFreq [Default no]
        -d <mindepth>          : minimum depth for BED and frCorr[Default 5]
===================================================================
Usage: bisReadMapper.pl -r <fq1[,fq2]> -W </path/to/WatsonIndex> -C </path/to/CrickIndex> -g </path/to/fai> -a <aligner> [options] > log
Required:
   -r   read(s) fq1[,fq2]
   -W   path to Watson converted index OR single combined W/C index
   -C   path to Crick converted index
   -g   path to indexed reference genome fasta file (.fai)
   -a   path to aligner, SOAP2 or Bowtie2 or path to LAST directory or path to GEM binary directory only
Optional settings:
   -e   encoded [1/0]
        If yes, then reads will not be encoded nor trimmed
   -p   number of processors for mapping [int], default 2
   -n   name of sample
   -3   3-prime trimming
   -b   ASCII base quality offset
   -5   5-prime trimming
   -q   quality value to perform quality trimming at
   -T   temporary directory for unix sort
   -s   --min-score function for bowtie2 (-L,-0.6,-0.6)
====================================================================

