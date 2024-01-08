#!/usr/bin/perl -w
use strict; use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd 'abs_path';
my ($marker,$sample,$outdir);
my (%clust,%gene);
GetOptions(
	"help|?" =>\&USAGE,
	"o:s"=>\$outdir,
	"marker:s"=>\$marker,
	"sample:s"=>\$sample,
) or &USAGE;
&USAGE unless ($outdir and $marker);

open MARKER,$marker or die $!;
<MARKER>;
while(<MARKER>){
	chomp;
	next if (/^,/);
	my ($id,$symbol,$pval,$logFC,$pct1,$pct2,$padj,$cluster,$gene,$regulate,$sig,$anno) = (split/,|\t/,$_,12)[0,1,2,3,4,5,6,7,8,9,10,11];
	#my ($id,$pval,$logFC,$pct1,$pct2,$padj,$cluster,$gene,$regulate,$sig,$anno) = (split/\t/,$_,11)[0,1,2,3,4,5,6,7,8,9,10];
	$clust{$cluster}{$gene} =  "$pct1\t$pct2" . "\t$logFC" . "\t$pval\t$padj\t$regulate\t$sig" . "\t$anno";
	$gene{$gene} = $id;
}
close MARKER;

my $annotation= "GO:biological_process\tGO:cellular_component\tGO:molecular_function\tKEGG:KO\tKEGG:Description";
system("rm -rf $outdir/$sample\_cluster_anno && mkdir -p $outdir/$sample\_cluster_anno");
for my $clus (sort keys %clust){
	open OUT,">","$outdir/$sample\_cluster_anno/$sample\_cluster$clus.anno.xls" or die $!;
	print OUT "GeneID\tSymbol\tCluster$clus\_pct\tCluster$clus\_Other_pct\tLogFoldChange\tPval\tAdjustPval\tUp/Down\tSignificant\t$annotation\n";
	for my $geneid(sort keys %{$clust{$clus}}){
		print OUT "$gene{$geneid}\t$geneid\t$clust{$clus}{$geneid}\n";
	}
	close OUT;
}

sub USAGE {
	my $usage=<<"USAGE";
#--------
Program: $0
Contact:junwuliu\@genome.cn
#--------
USAGE
	print $usage;
	exit 1;
}
