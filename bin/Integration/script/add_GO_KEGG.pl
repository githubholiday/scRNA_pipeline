#!/usr/bin/perl -w
use strict; use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd 'abs_path';
my ($in,$out,$anno,@title);
GetOptions(
	"help|?" =>\&USAGE,
	"out:s"=>\$out,
	"in:s"=>\$in,
	"anno:s" =>\$anno,
	"title=s{,}" =>\@title,
) or &USAGE;

&USAGE unless ($in and $out and $anno);

my %anno;
my ($add_title,$length) = &ReadAnno($anno);

open OUT,">",$out or die $!;
open IN,$in or die $!;
chomp(my $head = <IN>);
$head =~ s/,/\t/g;
print OUT "$head\t$add_title\n";
while(<IN>){
	chomp;
	my ($gene_id,$symbol,$description) = (split/,|\t/,$_,3)[0,1,2];
	$description =~ s/,/\t/g;
	my $annotation;
	if ($anno{$gene_id}){
		$annotation =  $anno{$gene_id};
	}elsif($anno{$symbol}){
		$annotation =  $anno{$symbol};
	}else{
		$annotation =  ".\t" x $length ;
		$annotation =~ s/\t$//;
	}
	#print($annotation) ; <STDIN>;
	print OUT "$gene_id\t$symbol\t$description\t$annotation\n";
}
close IN;

sub ReadAnno{
	my (@indexs,$add_title,$length);
	open ANNO,$anno or die $!;
	chomp(my $head = <ANNO>);
	my @header=split/\t/,$head;
	if (@title){
		for (@title){
			my $name = $_;
			my @index = grep {$header[$_] ~~ $name} 0..$#header;
			if (@index){
				push @indexs,$index[0];
			}else{
				die "$name is not found in $anno, please Check it\n";
			}
		}
		$add_title = join("\t",@title);
		$length = $#title + 1;
	}else{
		for (my $i=1;$i<=$#header;$i++){
			push @indexs,$i;
		}
		shift @header;
		my @add_title = @header;
		$add_title = join("\t",@add_title);
		$length = $#add_title + 1;
	}
	while(<ANNO>){
		chomp;
		my @desc = split/\t/;
		my $description;
		for (@indexs){
			$description .= $desc[$_] . "\t";
		}
		$description =~ s/\t$//;
		$anno{$desc[0]} = $description;
	}
	close ANNO;
	return ($add_title,$length);
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
