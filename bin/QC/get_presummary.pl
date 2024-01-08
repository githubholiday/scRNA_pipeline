#!/usr/bin/perl -w
use strict;

unless(@ARGV >=1){
	print "perl $0 <summary> <new_summary> <short_summary>\n";
	exit;
}


my$infile=$ARGV[0];
my$outfile=$ARGV[1];
my$shortfile=$ARGV[2];

open IN,$infile or die $!;
open OU,">$outfile" or die $!;
open SH,">$shortfile" or die $!;
my @short_list = ("Sample","Estimated_Number_of_Cells","Mean_Reads_per_Cell","Median_Genes_per_Cell","Number_of_Reads","Valid_Barcodes","Sequencing_Saturation","Reads_Mapped_to_Genome","Reads_Mapped_Confidently_to_Transcriptome","Fraction_Reads_in_Cells");
while(<IN>){
	chomp;
	$_=~s/%//g;
	my@tmp=split(/\t/,$_);
	$tmp[0]=~s/\s+/_/g;
	print OU join"\t",$tmp[0],@tmp[1..$#tmp];
	print OU "\n";
	if ( $tmp[0] ~~ @short_list ){
		print SH join"\t",$tmp[0],@tmp[1..$#tmp];
		print SH "\n";
	}
}
close IN;
close OU;
close SH;
