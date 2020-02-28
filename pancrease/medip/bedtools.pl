use strict;
my @sam=qw/2019032901 2019032903 2019040901 2019051703 2019052301 2019053101 2019053102/;
foreach my $sam(@sam){
my @fileC=glob("$sam\_N*narrowPeak");
my @fileT=glob("$sam\_T*narrowPeak");
foreach my $file(@fileT){
my ($out)=split/\_2019/,$file;
print "bedtools intersect -a $file -b $fileC[0] -v > $out.venn\n";
}
}
