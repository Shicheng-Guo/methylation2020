#!/usr/bin/perl

use strict;
use warnings;

use Algorithm::SVM::DataSet;
use Algorithm::SVM;

print("Creating new Algorithm::SVM\n");
my $svm = new Algorithm::SVM(Type => 'C-SVC',
				Kernel => 'radial');

print("Creating new Algorithm::SVM::DataSet objects\n");
my $ds1 = new Algorithm::SVM::DataSet(Label => 1);
my $ds2 = new Algorithm::SVM::DataSet(Label => 1);
my $ds3 = new Algorithm::SVM::DataSet(Label => 1);
my $ds4 = new Algorithm::SVM::DataSet(Label => 2);
my $ds5 = new Algorithm::SVM::DataSet(Label => 2);
my $ds6 = new Algorithm::SVM::DataSet(Label => 2);

print("Adding attributes to Algorithm::SVM::DataSet objects\n");
my @d1 = (1,0,1,1,1);
my @d2 = (1,1,0,1,1);
my @d3 = (1,1,1,0,1);
my @d4 = (0,1,0,0,0);
my @d5 = (1,0,0,0,0);
my @d6 = (1,0,1,1,1);

$ds1->attribute($_, $d1[$_ - 1]) for(1..scalar(@d1));
$ds2->attribute($_, $d2[$_ - 1]) for(1..scalar(@d2));
$ds3->attribute($_, $d3[$_ - 1]) for(1..scalar(@d3));
$ds4->attribute($_, $d4[$_ - 1]) for(1..scalar(@d4));
$ds5->attribute($_, $d5[$_ - 1]) for(1..scalar(@d5));
$ds6->attribute($_, $d6[$_ - 1]) for(1..scalar(@d6));

print("Training model\n");
my @tset=($ds1,$ds2,$ds3,$ds4,$ds5);
$svm->train(@tset);
my $v = $svm->predict_value($ds2);
my $data = $svm->save('sample.model.train2');

print("Checking predictions\n");

print("Here is accuracy: $v\n");
my $p1 = $svm->predict($ds1);
my $p2 = $svm->predict($ds2);
my $p3 = $svm->predict($ds3);
my $p4 = $svm->predict($ds4);
my $p5 = $svm->predict($ds5);
my $p6 = $svm->predict($ds6);

print("Here are the p1,p2,p3: $p1,$p2,$p3\n");
print("Here are the p4,p5,p6: $p4,$p5,$p6\n");

my $prob4 = $svm->getSVRProbability();
my $prob5 = $svm->getSVRProbability();
my $prob6 = $svm->getSVRProbability();



$svm=undef; # destroy svm object
