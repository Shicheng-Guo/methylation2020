use strict;
use Cwd;
open F, "filename.txt";
while(<F>){
my($srr,$filename)=split/\s|,/;
my $input="$srr.sorted.bam";
my $output="$filename.bam";
system("cp $input $output");
}
