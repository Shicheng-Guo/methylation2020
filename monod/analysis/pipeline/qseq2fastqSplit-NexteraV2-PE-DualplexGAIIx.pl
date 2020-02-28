#!/usr/bin/perl

use warnings;
use strict;
my $baseCallDir = $ARGV[0];
my $lane = $ARGV[1];
my $outSeqDir = $ARGV[2];

my @list_bases = ("A", "T", "G", "C", "N");

if(!$ARGV[2]){
	print "Converting and spliting compressed qseq files to fastq files based on Nextera's V2 barcodes.\n";
	print "Usage:\n";
	print "qseq2fastqSplit-NexteraV2-PE-DualplexGAIIx.pl  [basecallDir]  [lane]  [outSeqDir] [accepted 2mm N5] \n";	
	die();
}

# file _2_ is N7
# file _1_ is N5

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

my %N7barcode2idx = ("TAAGGCGA", "N701",
"CGTACTAG", "N702",
"AGGCAGAA", "N703",
"TCCTGAGC", "N704",
"GGACTCCT", "N705",
"TAGGCATG", "N706",
"CTCTCTAC", "N707",
"CAGAGAGG", "N708",
"GCTACGCT", "N709",
"CGAGGCTG", "N710",
"AAGAGGCA", "N711",
"GTAGAGGA", "N712");

my %N5barcode2idx = ("TAGATCGC", "N501",
"CTCTCTAT", "N502",
"TATCCTCT", "N503",
"AGAGTAGA", "N504",
"GTAAGGAG", "N505",
"ACTGCATA", "N506",
"AAGGAGTA", "N507",
"CTAAGCCT", "N508");

# Allow for 1 mm
my %mmN7barcode2idx;
my %mmN5barcode2idx;

# Allow these 2 mm codes for N5
my %mmN5barcode2idx_2;
if($ARGV[3]){
	open(INFILE, "$ARGV[3]") || die("Accepted N5 codes is not found\n");
	while(my $line = <INFILE>){
		chomp($line);
		my @f = split /\t/, $line;
		next if($f[2] =~ m/,/);
		$mmN5barcode2idx_2{$f[0]} = $f[2];
	}
}


my %outFileHandles;		

foreach my $code (keys(%N7barcode2idx)){
	# Allow for 1 mm
	my @bases = split "", $code;
	for(my $i = 0; $i < scalar(@bases); $i++){
		my @tmp = @bases;
		foreach my $nuc (@list_bases){
			$tmp[$i] = $nuc;
			my $mm = join "", @tmp;
			$mmN7barcode2idx{$mm} = $N7barcode2idx{$code};
		}
	}
}

foreach my $code (keys(%N5barcode2idx)){
        # Allow for 1 mm
        my @bases = split "", $code;
        for(my $i = 0; $i < scalar(@bases); $i++){
                my @tmp = @bases;
		foreach my $nuc (@list_bases){
	                $tmp[$i] = $nuc;
        	        my $mm = join "", @tmp;
                	$mmN5barcode2idx{$mm} = $N5barcode2idx{$code};
		}
        }
}

foreach my $codeN7 (keys(%N7barcode2idx)){
	foreach my $codeN5 (keys(%N5barcode2idx)){
		my $handle =  $N5barcode2idx{$codeN5} . "_" . $N7barcode2idx{$codeN7};
	        my $outFileName = $outSeqDir . "/" . $lane . "_1_" . $handle . ".txt";
        	open($outFileHandles{$handle}->{'1'}, ">$outFileName") || die("Error in opening file $outFileName!\n");
	        $outFileName = $outSeqDir . "/" . $lane . "_2_" . $handle . ".txt";
        	open($outFileHandles{$handle}->{'2'}, ">$outFileName") || die("Error in opening file $outFileName!\n");
	}
}

my ($total_reads, $quality_filtered_reads, $decoded_reads_pm_N7_N5, $decoded_reads_1mm_N7_N5, $decoded_reads_1mm_N7, $decoded_reads_1mm_N5, $decoded_reads_accepted_N5);
#my @tiles = (11,12,13,21,22,23);
#foreach my $i (@tiles) {	
for(my $j=1; $j<=120; $j++){
		my $seqFileName1 = sprintf("%s/%s_1_%04d_qseq.txt", $baseCallDir, $lane, $j);	
		open(SEQIN1, "$seqFileName1")||	die("Error in opening file $seqFileName1!\n");	
		my $seqFileName2 = sprintf("%s/%s_4_%04d_qseq.txt", $baseCallDir, $lane, $j);
		open(SEQIN2, "$seqFileName2")||	die("Error in opening file $seqFileName2!\n");	
		my $idxFileName1 = sprintf("%s/%s_2_%04d_qseq.txt", $baseCallDir, $lane, $j);
		open(IDXIN1, "$idxFileName1")||	die("Error in opening file $idxFileName1!\n");	
		 my $idxFileName2 = sprintf("%s/%s_3_%04d_qseq.txt", $baseCallDir, $lane, $j);
 		open(IDXIN2, "$idxFileName2")||   die("Error in opening file $idxFileName2!\n");

		print "Processing $seqFileName1 ...\n";			
		while (my $seq1 = <SEQIN1>) {
			my $seq2 = <SEQIN2>;
			my $index1 = <IDXIN1>;
			my $index2 = <IDXIN2>;
			chomp($seq1);
			chomp($seq2);
			chomp($index1);
			chomp($index2);
			$total_reads++;
			my @idx_parts = split(/\t/, $index1);
			next if($idx_parts[10] == 0);
			my $index_seq1 = $idx_parts[8];
			@idx_parts = split(/\t/, $index2);
			next if($idx_parts[10] == 0);
			my $index_seq2 = $idx_parts[8];
			$quality_filtered_reads++;
			my ($index_id1, $index_id2) = (0,0);
			my ($pm5, $pm7) = (0,0);
			if($N7barcode2idx{$index_seq1}){
				$index_id1 = $N7barcode2idx{$index_seq1};
				$pm7 = 1;
			}else{
				my @bases = split "", $index_seq1;
				for(my $i = 0; $i < scalar(@bases); $i++){
					my @tmp = @bases;
					$tmp[$i] = "N";
					my $mm = join "", @tmp;
					if($mmN7barcode2idx{$mm}){
						$index_id1 = $mmN7barcode2idx{$mm};
						last;
					}
				}
				next if(!$index_id1);
			}
			if($N5barcode2idx{$index_seq2}){
                                $index_id2 = $N5barcode2idx{$index_seq2};
                                $pm5 = 1;
                        }else{
                                my @bases = split "", $index_seq2;
                                for(my $i = 0; $i < scalar(@bases); $i++){
                                        my @tmp = @bases;
                                        $tmp[$i] = "N";
                                        my $mm = join "", @tmp;
                                        if($mmN5barcode2idx{$mm}){
                                                $index_id2 = $mmN5barcode2idx{$mm};
                                                last;
                                        }
                                }
                                if(!$index_id2){
					#use 2 mm
					if($mmN5barcode2idx_2{$index_seq2}){
						$index_id2 = $mmN5barcode2idx_2{$index_seq2};
						$pm5 = -1;
					}else{
						next;
					}
				}
                        }

			$decoded_reads_pm_N7_N5++ if($pm5 == 1 && $pm7 == 1);
			$decoded_reads_1mm_N7_N5++ if($pm5 == 0 && $pm7 == 0);
			$decoded_reads_1mm_N7++ if($pm5 == 1 && $pm7 == 0);
			$decoded_reads_1mm_N5++ if($pm5 == 0 && $pm7 == 1);
			$decoded_reads_accepted_N5++ if($pm5 == -1);
			my $handle = $index_id2 . "_" . $index_id1;
			my @parts = split(/\t/, $seq1);
			$parts[8] =~ s/\./N/g;
			my $fileHandle = $outFileHandles{$handle}->{'1'};
			print $fileHandle "@","$parts[0]:$parts[2]:$parts[3]:$parts[4]:$parts[5]#$parts[6]/$parts[7]\n";
			print $fileHandle "$parts[8]\n";
			print $fileHandle "+","$parts[0]:$parts[2]:$parts[3]:$parts[4]:$parts[5]#$parts[6]/$parts[7]\n";
			print $fileHandle "$parts[9]\n";
			
			@parts = split(/\t/, $seq2);
			$parts[8] =~ s/\./N/g;
			$fileHandle = $outFileHandles{$handle}->{'2'};
			print $fileHandle "@","$parts[0]:$parts[2]:$parts[3]:$parts[4]:$parts[5]#$parts[6]/$parts[7]\n";
			print $fileHandle "$parts[8]\n";
			print $fileHandle "+","$parts[0]:$parts[2]:$parts[3]:$parts[4]:$parts[5]#$parts[6]/$parts[7]\n";
			print $fileHandle "$parts[9]\n";			
		}
		print "Total: $total_reads\t",
			"Passed QC: $quality_filtered_reads\t",
			"Decoded(pm): $decoded_reads_pm_N7_N5\t",
			"Decoded(1mm N7): $decoded_reads_1mm_N7\t",
			"Decoded(1mm N5): $decoded_reads_1mm_N5\t",
			"Decoded(accepted N5): $decoded_reads_accepted_N5\t",
			"Decoded(1mm both): $decoded_reads_1mm_N7_N5\n";
	close(SEQIN1);
	close(SEQIN2);
	close(IDXIN1);
	close(IDXIN2);
}

foreach my $index_id (keys(%outFileHandles)){
	close($outFileHandles{$index_id}->{'1'});	
	close($outFileHandles{$index_id}->{'2'});	
}

sub revComp(){
	my $seq = shift;
	my $rcSeq='';
	for(my $i=0; $i<length($seq); $i++){
		$rcSeq = $rcTable{substr($seq,$i,1)} . $rcSeq;
	}
	return $rcSeq;
}
