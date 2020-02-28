# This software is Copyright © 2017 The Regents of the University of California. All Rights Reserved.
#  
# Permission to copy, modify, and distribute this software and its documentation for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
#  
# Permission to make commercial use of this software may be obtained by contacting:
# Office of Innovation and Commercialization
# 9500 Gilman Drive, Mail Code 0910
# University of California
# La Jolla, CA 92093-0910
# (858) 534-5815
# kzhang@ucsd.edu
 
# This software program and documentation are copyrighted by The Regents of the University of California. The software program and documentation are supplied "as is", without any accompanying services from The Regents. The Regents does not warrant that the operation of the program will be uninterrupted or error-free. The end-user understands that the program was developed for research purposes and is advised not to rely exclusively on the program for any reason.
#  
# IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
# CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
# EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF
# CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
# MODIFICATIONS.

# this script reads in a haploInfo file and bed regions and it invoke hapinfo2LDR2.pl
# Output the average rsq value within the bed regions
# Run the script to the Bam directory
# Contact: Shicheng Guo <shihcheng.guo@gmail.com>
# Version 1.3
# Update: 2016-02-25

#!/usr/bin/perl
use Cwd;
chdir getcwd;
printUsage() if scalar(@ARGV)<2;
my $bed=shift @ARGV;
my $hapinfo=shift @ARGV;
open F1,$bed || die "cannot open $bed\n";
while(<F1>){
chomp;
next if /^\s+$/;
my ($chr,$start,$end)=split/\t/;
my $cor="$chr:$start-$end";
system("perl /home/shg047/bin/hapinfo2LDR2.pl tmp $cor \< $hapinfo");
open F2,"tmp.$chr.rsq" || die "cannot open tmp.$chr.rsq\n";
my $i=0;
my $Rlt=0;
my $output;
# print "----------------------------------------------------------\n";
while(<F2>){
# print "$_";
next if /^\t/;
my (undef,@line)=split/\s+/;
foreach my $tmp(@line){
next if $tmp=~/NA/i;
#print("$i\t$tmp\t$Rlt\n");
$i++;
$Rlt=$Rlt+$tmp;	
}
}
# print "-----------------------------------------------------------\n";
if($i gt 1){
$output=$Rlt/($i);
}else{
$output="NA";
}

print "$cor\t$output\n";
}


sub printUsage{
        print " \nUsage:\n";
                print " perl $0 input.bed input.hapinfo\n";
                print " For example:\n perl $0 chr1.bed haploinfo.input.txt\n\n";
                print "-----------------------------------------------------------------------------\n";
                print " haploinfo.input.txt format:\n";
                print " chr10:10000873-10001472        CCC     1       10001056,10001082,10001168\n\n";
        print "-----------------------------------------------------------------------------\n";
        print " Output Format:\n";
        print " chr10:10000873-10001472\t0.64\n\n";
        print "-----------------------------------------------------------------------------\n";

        exit 0;
}
