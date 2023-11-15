#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my ($result,$pval,$qval,$lgfc,$out,$tsv);
my (%gene);
GetOptions(
	"h|help|?" => \&USAGE, 
	"result:s" => \$result,
	"p:s" => \$pval,
	"q:s" => \$qval,
	"lgfc:s" => \$lgfc,
	"out:s" => \$out,
	"tsv:s" => \$tsv,
);

$pval ||= "0.05";
$qval ||= "0.05";
$lgfc ||= "0.25";

&ReadTSV;

open OUT,">",$out or die $!;
open FL,$result or die $!;
my $title = <FL>;
chomp($title);
$title =~ s/^,/Gene_Symbol,/;
#$title =~ s/^,//;
print OUT "Gene_ID,$title,Up/Down,Significant\n";
while(<FL>){
	chomp;
	my ($sig,$regulate);
	my $line = $_;
	my ($id,$p,$logfc,$pct1,$pct2,$q,$cluster) = (split/,|\t/)[0,1,2,3,4,5,6];
	#print "$logfc\t$p\t$q\n";
	$q = "1" if ($q eq "NA");
	$p = "1" if ($p eq "NA");
	$regulate = ($logfc > 0) ? "Up" : "Down";
	if ( (abs($logfc) >= $lgfc) && $p <= $pval && $q <= $qval ){
		$sig = "yes";
	}else{
		$sig = "no";
	}
	print STDERR if (!$gene{$id});
	#$line =~s/^$gene,// ;
	print OUT "$gene{$id},$line,$regulate,$sig\n";
}
close FL;
close OUT;

sub ReadTSV{
	open TSV,"gzip -dc $tsv |" or die $!;
	while(<TSV>){
		chomp;
		next if /^$/;
		my ($id,$name) = (split/\t/)[0,1];
		$gene{$name} = $id;
	}
	close TSV;
}

sub USAGE{
	my $usage=<<"USAGE";
----------------------
Contact: junwuliu\@annoroad.com
   Date: 2016-12-18
  Usage:
	perl $0 -i sample1_sample2.report.xls -p 0.05 -q 0.05 -fc 2 -o sample1_sample2.report_adjust.xls

----------------------
USAGE
	print $usage;
	exit;
}
