#!/usr/bin/perl -w

use strict;

my $samplesheet = $ARGV[0];
my $number_mm = $ARGV[1];
my $r1_file = $ARGV[2];
my $r2_file = $ARGV[3];

my %IndxTable;
my %fileHandlers;

my $total_reads;
my $decoded_reads;
my $decoded_0mm_reads = 0;
my $decoded_1mm_reads = 0;
my $decoded_2mm_reads = 0;
my $decoded_3mm_reads = 0;
open(IN, "$samplesheet") || die("Error opening file\n");
my $line = <IN>;
$line = <IN>;
while($line = <IN>){
	#Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description

	chomp($line);
	my @f = split /,/, $line;
	my ($id, $barcode7, $barcode5) = ($f[1], $f[6],$f[8] );
	print "Original: $id \t", $barcode7, $barcode5,"\n";
	$IndxTable{$barcode7.$barcode5} = $id;
	# make fileHandlers
	open( $fileHandlers{$id}->{"r1"} , ">$id.R1.fastq" ) || die("Error oepning $id.R1.fastq");
	open( $fileHandlers{$id}->{"r2"} , ">$id.R2.fastq" ) || die("Error oepning $id.R2.fastq") if($r2_file);
}
open(IN1, "$r1_file") || die("Error opening $r1_file\n");
open(IN2, "$r2_file") || die("Error opening $r2_file\n") if($r2_file);
while(my $line1 = <IN1>){
	my $line2 = <IN1>;
	my $line3 = <IN1>;
	my $line4 = <IN1>;
	my $line1_2 = <IN2>;
	my $line2_2 = <IN2>;
	my $line3_2 = <IN2>;
	my $line4_2 = <IN2>; 
	my ($id, $unknownbarcode) = split " ", $line1;
	$unknownbarcode =~ s/1:N:0://g;
	$unknownbarcode =~ s/[+]//g;
	#print $unknownbarcode, "\n";
	my $min_edit_distance = 1000;
	my @ids_with_min_distance;
	my @unk = split "", $unknownbarcode;
	foreach my $barcode(keys %IndxTable){
		my $cur_diff = 0;
		my @bases = split "", $barcode;
		for(my $j = 0; $j < scalar(@bases); $j++){
			if($unk[$j] ne $bases[$j]){
				$cur_diff++;
			}
		}
		if($cur_diff < $min_edit_distance){
			$min_edit_distance = $cur_diff;
			while(scalar(@ids_with_min_distance) > 0){
				pop(@ids_with_min_distance);
			}
			push(@ids_with_min_distance, $IndxTable{$barcode});
		}elsif($cur_diff eq $min_edit_distance){
			push(@ids_with_min_distance, $IndxTable{$barcode});
		}
	}
	$total_reads++;
	next if($min_edit_distance > $number_mm);
	if(scalar(@ids_with_min_distance) == 1){
		$decoded_0mm_reads++ if($min_edit_distance==0);
		$decoded_1mm_reads++ if($min_edit_distance==1);
		$decoded_2mm_reads++ if($min_edit_distance==2);
		$decoded_3mm_reads++ if($min_edit_distance==3);
		$decoded_reads++;
		print {$fileHandlers{$ids_with_min_distance[0]}->{"r1"}} $line1,$line2,$line3,$line4;
		print {$fileHandlers{$ids_with_min_distance[0]}->{"r2"}} $line1_2,$line2_2,$line3_2,$line4_2;
	}
}
#close(IN1);
#close(IN2);

foreach my $id (keys %fileHandlers){
	close( $fileHandlers{$id}->{"r1"} );
	close( $fileHandlers{$id}->{"r2"} );
}

print "Total reads: $total_reads\n";
print "Decoded reads: $decoded_reads/$total_reads ", sprintf("%3.2f", $decoded_reads/$total_reads), "\n";
print "0 mm reads: $decoded_0mm_reads/$decoded_reads ", sprintf("%3.2f", $decoded_0mm_reads/$decoded_reads), "\n";
print "1 mm reads: $decoded_1mm_reads/$decoded_reads ", sprintf("%3.2f", $decoded_1mm_reads/$decoded_reads), "\n";
print "2 mm reads: $decoded_2mm_reads/$decoded_reads ", sprintf("%3.2f", $decoded_2mm_reads/$decoded_reads), "\n";
print "3 mm reads: $decoded_3mm_reads/$decoded_reads ", sprintf("%3.2f", $decoded_3mm_reads/$decoded_reads), "\n";
