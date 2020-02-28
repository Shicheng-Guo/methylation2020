#!/usr/bin/perl -w
use strict;

# INPUT
# lane#	project	i5_file i7_file id_file

# OUTPUT
# Lane Sample_ID	Sample_Name	Sample_Plate	Sample_Well	I7_Index_ID	index	I5_Index_ID	index2	Sample_Project	Description
print "[Data]\n";
print "Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description\n";

open(TABLEFILE, "$ARGV[0]") || die("Error reading input table");
while(my $fline = <TABLEFILE>){
	chomp($fline);
	my @tmp = split "\t", $fline;
	my %i7_values;
	my %i5_values;
	my %ids;
	if($tmp[4]){
		my ($lane, $project) = ($tmp[0], $tmp[1]);
		my ($i5_file, $i7_file, $id_file) = ($tmp[2], $tmp[3], $tmp[4]);
		open(INFILE, "$id_file") || die("Error opening $id_file\n");
		while(my $line = <INFILE>){
	        	chomp($line);
		        my ($sid, $sname) = split "\t", $line;
        		$ids{$sname} = $sid;
		}
		close(INFILE);
		open(IN5, "$i5_file") || die("Error opening $i5_file\n");
		open(IN7, "$i7_file") || die("Error opening $i7_file\n");
		while(my $line = <IN5>){
			chomp($line);
			my ($id, $seq) = split "\t", $line;
			$i5_values{$id} = $seq;
		} 
		while(my $line = <IN7>){
			chomp($line);
			my ($id, $seq) = split "\t", $line;
			$i7_values{$id} = $seq;
		}
		foreach my $i7 (keys %i7_values){
			foreach my $i5 (keys %i5_values){
				my $sname = "$i7-$i5";
				if($ids{$sname}){
					my $sid = $ids{$sname};
					$sname = $sid;
				}
				my ($seq7, $seq5) = ($i7_values{$i7}, $i5_values{$i5});
				print "$lane,$sname,,,,$i7,$seq7,$i5,$seq5,$project,\n";
			}
# Lane Sample_ID	Sample_Name	Sample_Plate	Sample_Well	I7_Index_ID	index	I5_Index_ID	index2	Sample_Project	Description
		}
	}
	else{
		my $in_file = $tmp[2];
		my $lane = $tmp[0];
		my $project =$tmp[1];
		open(IN, "$in_file") || die("Error opening $in_file\n");
		while(my $line = <IN>){
			chomp($line);
			my ($id, $seq) = split "\t", $line;
			print "$lane,$id,,,,,$seq,,,$project,\n";
		}
	}
	foreach my $value(keys %i7_values){
		delete($i7_values{$value});
	}
	foreach my $value(keys %i5_values){
		delete($i5_values{$value});
	}
	foreach my $value(keys %ids){
		delete($ids{$value});
	}
}
close(TABLEFILE);
