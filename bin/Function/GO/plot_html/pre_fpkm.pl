#!/usr/bin/perl -w
use strict;

unless(@ARGV>=1){
	print "perl $0 <fpkm.xls> <header> <fpkm.tmp.xls>\n";
	exit;
}

my $infile=$ARGV[0];
my $header=$ARGV[1];
my $outfile=$ARGV[2];

open IN,$infile or die $!;
open OUT,">$outfile" or die $!;

if($header eq "T"){
	my$title = <IN> ;
	chomp $title;
	print OUT "$title\n";
}
else{ print "No header,go on!";}	
while(<IN>){
	chomp;
	my@matrix=split(/\t/,$_);
	print OUT"$matrix[0]";
	for(my$i = 1;$i < @matrix ;$i++){
		my$tmp=$matrix[$i]+1;
		my$log=log($tmp)/log(10);
		print OUT"\t$log";	
	}
	print OUT"\n";	
}
close IN;
