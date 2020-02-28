use strict;
my @file=glob("*.o*");
my %data;
foreach my $file(@file){
if($file =~ /(\d+_\w\d*)_2019.+\.(\d+_\w\d*)_2019./i){
open F,$file;
while(<F>){
my ($r2)=split/\s+/;
$data{$1}{$2}=$r2;
}
}
}

my @s=sort keys %data;
my $header=join("\t",@s);
print "\t$header\n";
foreach my $s1(sort keys %data){
   print "$s1";
   foreach my $s2(sort keys %data){
   print "\t$data{$s1}{$s2}";
   }
   print "\n";
}

