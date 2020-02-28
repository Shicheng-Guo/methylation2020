#!/usr/bin/perl 
use strict;
use warnings;
my $dimension;my $data;   #A reference to three dimesion array
#Initialization #3x3x3 #Value range from 000 to 333
my $dir="/work/w¹¤×÷»ã±¨/meta/datasets/TCGALung";
chdir $dir or die"";
###############################################to get gene list
my %batchs;
my $genefile="/work/data/methylome/TCGA/zMeth27Prob/probehg19.bed"; #Methylation Gene list
open F,"$genefile";
my @gene;
while (<F>){
chomp;
next if /gene/;
my($cg)=split/\t/,$_;
push @gene,$cg;
}
close F;
##############################################to get individual list
open OUT,">Methylation.txt";
print OUT "\t";

system("ls *TCGA* >IndividualList.txt");
open F,"IndividualList.txt";
while (<F>){
 my $line=$_;
 $line=~/(TCGA-\w+(-\d+)+)|(TCGA(-\d+)+)/g;  
 my $sb=$1;
 print "$sb\n";
 $sblist{$sb}=$sb;	
}
close F;
##############################################to intergret Methylation or expression Datasets
print "Now Reading array information,Please Waiting \n";
my @files=glob "*TCGA*.txt";
for my $files(@files){
	$files=~/(TCGA-\w+(-\d+)+)|(TCGA(-\d+)+)/g;
    my $indi=$1;
	open F,$files;
	my $i;	 
	while(<F>){
	    chomp; chop;	   
            next if /gene|barcode|Composite|REF/;
	    $i++;
	    my @line=split /\t/;			
	    if ($files=~m/methylation/i){
	    	$data{$indi}{$line[0]} = $line[1];
	    }	
}
}


print "Now writing array information to OutputFiles,Please Waiting\n";
my($j,$i,$k);
for my $cg(@cg){
print OUT "$cg\t";
}
print OUT "\n";


foreach (keys %sblist){
   $i = $_;
   my (undef,undef,undef,$sampletype)=split /-/,$i;
   if(!defined $batchs{$i}){ print "$i Batch Information Lost,Please Check\n";}
   print OUT "$i\t$sampletype\t$sblist2{$i}\t$batchs{$i}\t";
   print OUT2 "$i\t";
   print OUT3 "$i\t";
   
   foreach my $gene( @gene ){
      ($j,undef) = split(/_/,$gene);
      foreach(@type){
      	 $k=$_;
      	if($k eq "methylation"){
          if(defined $data{$i}{$gene}{$k} && $data{$i}{$gene}{$k} ne "null" && $data{$i}{$gene}{$k} ne "NULL"){
          print OUT "$data{$i}{$gene}{$k}\t";
          print OUT2 "$data{$i}{$gene}{$k}\t";
      #    print  "$data{$i}{$gene}{$k}\t";

          }else{
          print OUT "NA\t";
          print OUT2 "NA\t";
         }     
      	}else{
      		if (defined $data{$i}{$j}{$k} && $data{$i}{$j}{$k} ne "null" && $data{$i}{$j}{$k} ne "NULL"){
      			print OUT "$data{$i}{$j}{$k}\t";
      			print OUT3 "$data{$i}{$j}{$k}\t";
      	#		print  "$data{$i}{$j}{$k}\t";
      			
      			
      		}else{
      			print OUT "NA\t";
                print OUT3 "NA\t";
      		}
     
      	}
      }
   }
     print OUT "\n" ;
     print OUT2 "\n" ;
     print OUT3 "\n" ;
}

print "OK!\n";
print "The Output is MeEXP.txt";
close OUT;
close OUT2;
close OUT3;
