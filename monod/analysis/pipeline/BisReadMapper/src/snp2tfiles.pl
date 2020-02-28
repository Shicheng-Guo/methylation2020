#!/usr/bin/perl -w

use strict;

my %three2one;
$three2one{"A/G"}="R";
$three2one{"C/T"}="Y";
$three2one{"A/C"}="M";
$three2one{"G/T"}="K";
$three2one{"C/G"}="S";
$three2one{"A/T"}="W";

my %one2two;
$one2two{"A"} = "A A";
$one2two{"T"} = "T T";
$one2two{"G"} = "G G";
$one2two{"C"} = "C C";
$one2two{"R"} = "A G";
$one2two{"Y"} = "C T";
$one2two{"M"} = "A C";
$one2two{"K"} = "G T";
$one2two{"S"} = "C G";
$one2two{"W"} = "A T";

my %snp;
my %rs;
my $name = $ARGV[0];
while(my $line = <STDIN>){
	next if($line =~ m/SNP/);
	my @f = split /\t/, $line;
	#chr16:85185427  G       108     rs924475        S       G       G,27    -      -
	next if($f[1] eq "?");
	my $call = $one2two{$f[1]};
	if($f[3] ne "-"){
		$rs{$f[0]} = $f[3];
	}
	if($snp{$f[0]}){
		$snp{$f[0]} = "bad" if($call ne $snp{$f[0]});
	}else{
		$snp{$f[0]} = $call;
	}
}

my $id = $name;
$id =~ s/Indx//;
open(FAM, ">$name.tfam") || die("Error opening fam file\n");
print FAM "$id $name 0 0 0 -9\n";
close(FAM);
open(PED, ">$name.tped") || die("Error opening ped file\n");
foreach my $key(keys %snp){
	next if($snp{$key} eq "bad");
	my ($chr, $pos) = split ":", $key;
	my $rs_value = $chr.":".$pos;
	$rs_value = $rs{$key} if($rs{$key});
	my $num_c = $chr;
	$num_c =~ s/chr//g;
	print PED "$num_c $rs_value 0 $pos ", $snp{$key}, "\n";
}
close(PED);
