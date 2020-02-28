#!/usr/bin/perl -w
my $fq1 = $ARGV[0];
my $fq2 = $ARGV[1];
my $name = $ARGV[2];
my $umi_length = $ARGV[3];

my $se=0;

if($fq1 =~ /\.gz$/){
	open(IN1, "gunzip -c $fq1 |") || die("Error opening $fq1\n");
}else{
	open(IN1, "$fq1") || die("Error opening $fq1\n");
}

if($ARGV[3]){
	if($fq2 =~ /\.gz$/){
		open(IN2, "gunzip -c $fq2 |") || die("Error opening $fq2\n");	
	}else{
		open(IN2, "$fq2") || die("Error opening $fq2\n");
	}
}else{$se = 1; $umi_length = $ARGV[2]; $name = $ARGV[1];}
open(OUT1, ">$name.R1.fq") || die("Error opening $name.R1.fq\n");
open(OUT2, ">$name.R2.fq") || die("Error opening $name.R2.fq\n") if($se eq 0);

while(my $line1 = <IN1>){
	my $line2 = <IN1>;
	my $line3 = <IN1>;
	my $line4 = <IN1>;
	chomp($line1);
	chomp($line2);
	chomp($line3);
	chomp($line4);
	$line1 =~ s/^.//s;
	my @f = split /[\s+\t]/, $line1;
	$line1 = join("_", @f);
	my $umi = substr($line2, 0, $umi_length);
	my $s = substr($line2, ($umi_length), length($line2)-($umi_length));
	my $q = substr($line4, ($umi_length), length($line4)-($umi_length));
	print OUT1 "\@",$umi, "_" , "$line1\n$s\n+\n$q\n";
	next if($se eq 1);
	my $line1b = <IN2>;
	my $line2b = <IN2>;
	my $line3b = <IN2>;
	my $line4b = <IN2>;
	chomp($line1b);
	chomp($line2b);
	chomp($line3b);
	chomp($line4b);
	$line1b =~ s/^.//s;
	@f = split /[\s+\t]/, $line1b;
	$line1b = join("_", @f);
	$s = $line2b;
	$q = $line4b;
	print OUT2 "\@",$umi, "_", "$line1b\n$s\n+\n$q\n";
}
close(IN1);
close(IN2) if($se eq 0);
close(OUT2) if($se eq 0);
close(OUT1);
