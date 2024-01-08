#!/usr/bin/perl -w
use strict; use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd 'abs_path';
my ($up,$down,$ko2map,$keg,$out,$ko);
GetOptions(
	"help|?" =>\&USAGE,
	"up:s"=>\$up,
	"down:s"=>\$down,
	"ko2map:s" => \$ko2map, ### ko2map.xls
	"keg:s" => \$keg, ### ko00001.keg
	"ko:s" => \$ko, ## ko.list, KEGG_annotate
	"out:s" => \$out,
) or &USAGE;
&USAGE unless ($up and $down and $ko2map and $out and $ko and $keg);
my (%keg,%ko,%map);

&ReadKeg;
&ReadKo2Map;
&ReadKO;

open OUT,">",$out or die $!;
print OUT "#Gene\tKO\tMap\tFeature\tDescription\n";
&ReadDE($up,"up");
&ReadDE($down,"down");
close OUT;

sub ReadDE{
	my ($file,$type) = @_;
	open FL,$file or die $!;
	while(<FL>){
		chomp;
		my ($gene) = (split/\t/)[0];
		if ($ko{$gene}){
			for my $K(sort keys %{$ko{$gene}}){
				my $map_id = $map{$K} ? $map{$K} : ".";
				my $map_des = $keg{$K} ? $keg{$K} : ".";
				print OUT "$gene\t$K\t$map_id\t$type\t$map_des\n";
			}
		}else{
			print OUT "$gene\t.\t.\t$type\t.\n";
		}
	}
	close FL;
}

sub ReadKO{
	open KO,$ko or die $!;
	while(<KO>){
		chomp;
		my ($gene,$k_id) = (split/\t/)[0,1];
		next if ($k_id eq ".");
		$ko{$gene}{$k_id} = 1;
	}
	close KO;
}

sub ReadKo2Map{
	open KO,$ko2map or die $!;
	while(<KO>){
		chomp;
		my ($k_id,$map) = (split/\t/)[0,1];
		$map =~ s/\|/,/g;
		$map{$k_id} = $map;
	}
	close KO;
}

sub ReadKeg{
	open KEG,$keg or die $!;
	while(<KEG>){
		chomp;
		if (/(K\d+)\s+(.*)/){
			$keg{$1} = $2;
		}
	}
	close KEG;
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
