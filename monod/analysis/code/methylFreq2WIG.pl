#!/usr/bin/perl -w

use strict;

my $sampleID = $ARGV[0];
$sampleID = "Sample" if (!$sampleID);
my $minDepth = $ARGV[1];
$minDepth = 10 if(!$minDepth);

print "track type=wiggle_0 name=\"$sampleID:methRatio\" visibility=full color=20,150,20 altColor=150,20,20 windowingFunction=mean\n";

my %methylTable;
sub main(){
        while(my $line = <STDIN>){
                chomp($line);
                my @fields = split(/\t/, $line);
                my $strand = $fields[2] eq 'W' ? '+' : '-';
                my %alleleCounts;
                my $CT_counts;
                for(my $i=5; $i<scalar(@fields); $i+=2){
                        $alleleCounts{$fields[$i]}=$fields[$i+1];
                        $CT_counts += $fields[$i+1] if($fields[$i]=~ /[CT]/);
                }
                next if(!$CT_counts || $CT_counts/$fields[3] < 0.9);
                my $index=$fields[0] . ":" . $fields[1];
                $alleleCounts{'C'} =0 if(!$alleleCounts{'C'});
                $methylTable{$fields[0]}->{$fields[1]}->{'C'} +=  $alleleCounts{'C'} ;
                $methylTable{$fields[0]}->{$fields[1]}->{'CT'} += $CT_counts;
        }
        report_methylFreqBED();
}

sub report_methylFreqBED(){
        my $cur_chr = "NA";
        foreach my $chr(sort keys(%methylTable)){
		print "variableStep chrom=$chr\n";
		foreach my $pos(sort {$a<=>$b} keys %{$methylTable{$chr}}){
	                next if($methylTable{$chr}->{$pos}->{'CT'}<$minDepth);
                	my $methylLevel = sprintf("%4.3f", $methylTable{$chr}->{$pos}->{'C'}/$methylTable{$chr}->{$pos}->{'CT'});
			print $pos-1, "\t", $methylLevel, "\n";	
			print $pos, "\t", $methylLevel, "\n";
		}
        }
}

main();

