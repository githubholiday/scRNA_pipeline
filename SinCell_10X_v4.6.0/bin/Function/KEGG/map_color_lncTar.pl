#!/usr/bin/perl -w
use strict;

use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);
use Cwd qw(abs_path realpath getcwd);
use Getopt::Long;

my $version = "1.0";

###############################################################
#
# 	               set parameters,options
#
###############################################################
my ($kegg_report,$outdir,$mapdb); 

GetOptions(
	  'h|help|?' => \&USAGE,
	  'kegg_report:s' => \$kegg_report,
	  'mapdb:s' => \$mapdb,
	  'od:s' => \$outdir,

);

&USAGE unless ( $kegg_report && $mapdb ) ;

### softwares
my $wkhtmltoimage="$Bin/wkhtmltoimage-amd64-master/bin/wkhtmltoimage-amd64";

### default parameters
$outdir ||= "kegg_result";
$outdir = abs_path($outdir);

### 
if ( -d $outdir ) {
	system "rm -rf $outdir && mkdir -p $outdir";
}else{
	system "mkdir -p $outdir";
}

##############################################################
#
#		            kegg report file
#
##############################################################			
my %deg_info;
my %coords_info;
my %count;
my @filter_map = ("map01100","map01110","map01120","map01130");# these maps is very big and complicated

open KEGG_REPORT,$kegg_report or die "cannot find $kegg_report";
<KEGG_REPORT>;
while(<KEGG_REPORT>){
	chomp;
	my (@up_K_numbers,@down_K_numbers,@normal_K_numbers,@K_numbers,%count_of_up_down_gene_in_K);

	my ($mapid,$gene_de,$up_genes,$up_gene_count,$down_genes,$down_gene_count,$normal_genes) = (split /\t/, $_)[1,9,10,11,12,13,16];
	$normal_genes ||= "";
	#next if($up_gene_count+$down_gene_count == 0 );
	next if ($gene_de eq ".");
	my $map_html = "$mapdb/$mapid.html";
	my $map_png = "$mapdb/$mapid.png";
	
	# check if the map***.png exists
	unless ( -e $map_png ) {
		print " The $map_png not exists!\n";
		next;
	}

	# filter maps which is very big and complicated
	if ( grep { $_ eq $mapid } @filter_map ) {
		print "The $mapid is filter because it is a whole pathway map\n";
		next;
	}

	# kegg map containing up-regulated genes  
	unless ( $up_genes eq "." ) {	
		my @up_gene_K_numbers = split /;/,$up_genes;
		@up_K_numbers = map { $_."|up" } @up_gene_K_numbers;
	}

	# kegg map containing down-regulated genes
	unless ( $down_genes eq "." ) {
		my @down_gene_K_numbers = split /;/,$down_genes;
		@down_K_numbers =  map { $_."|down" } @down_gene_K_numbers;
	}
	unless ($normal_genes eq ""){
		my @normal_gene_K_numbers = split /;/,$normal_genes;
		@normal_K_numbers = map {$_."|normal"} @normal_gene_K_numbers;
	}

	#@K_numbers = (@up_K_numbers,@down_K_numbers);
	@K_numbers = (@up_K_numbers,@down_K_numbers,@normal_K_numbers);
	# test point 1
	#print "#merged K number\n";
	#print join("#",@K_numbers),"\n";

	# K_number is up-regulation or down-regulation or both ? 
	foreach my $K_number ( @K_numbers ) {
		my @arr = split /\|/,$K_number;
		foreach my $index ( 1..($#arr-1)){
			push @{$deg_info{$arr[$index]}},$arr[-1]; 
			$count_of_up_down_gene_in_K{$arr[$index]}->{$arr[-1]}++; 
		}
	}


	&rmdup(\%deg_info);
	
	# test point 2
	# print Dumper (\%deg_info);

	foreach my $key ( keys %count_of_up_down_gene_in_K ) { # judge K_numer is up-regulation or down-regulation or both by counts of gene up/down-regulation

		my @up_down_type = keys %{$count_of_up_down_gene_in_K{$key}};
		
		if ( @up_down_type >1 ) {
			my $up_count = $count_of_up_down_gene_in_K{$key}->{up} ? $count_of_up_down_gene_in_K{$key}->{up} : 0 ;
			my $down_count = $count_of_up_down_gene_in_K{$key}->{down} ? $count_of_up_down_gene_in_K{$key}->{down} : 0;
			my $normal_count = $count_of_up_down_gene_in_K{$key}->{normal} ? $count_of_up_down_gene_in_K{$key}->{normal} : 0;
			if ( $up_count > $down_count && $up_count > $normal_count) {
				@{$deg_info{$key}} = ("up");
			
			}elsif( $up_count < $down_count && $normal_count < $down_count ) {
				@{$deg_info{$key}} = ("down");
			
			}elsif ( $up_count == $down_count && $up_count > 0 ) {
				@{$deg_info{$key}} = ("down","up");
			}else{
				@{$deg_info{$key}} = ("normal");	
			}
		}
	}
	
	# test point 3
	#print Dumper(\%count_of_up_down_gene_in_K);
	
	
	### Get Gene Coords		
	&get_coords($mapid,$map_html,\%coords_info,\%count_of_up_down_gene_in_K);
	undef %count_of_up_down_gene_in_K;

	# test point 4
	#print Dumper(\%coords_info);

	### Color Map 
	&color_map_html($mapid,$map_png,$outdir);

	### end 
	undef %deg_info;
	undef %coords_info;
	print "The $mapid.png has been colored successfully!\n"; 

}
close KEGG_REPORT;

##############################################################
#
#		   	sub  functions
#
############################################################## 

sub rmdup{
	my ($hash) = @_;
	
	# remove duplicated value
	foreach my $key ( sort keys %{$hash} ) {
		my %count;
		my @uniq = sort grep { ++$count{$_} < 2 } @{$hash->{$key}};
		$hash->{$key} = \@uniq;
	}
}		

sub get_coords{
	my ($mapid,$map_html, $coords_info,$count_of_up_down_gene_in_K ) = @_;

	# load mapid.html
	open HTML,$map_html or die "cannot find $map_html";
	while(<HTML>){
		chomp;
		if ( $_ =~ /^<area shape=rect\tcoords=(\d+,\d+,\d+,\d+).*?title=\"(K\d+[^"]+)\"/ ){
			my $coords = $1;
			my $title = $2;
			my (@k_nums) = $title =~ /(K\d+) \((.+?)\)/g;
 			my ($enzyme) = $title =~ /, (\d+\.\d+\.\d+\.[\d-]+)/;


			if ( @k_nums == 2 ){# title is only K-number ( gene1 ), $enzyme,such as title="K01008 (selD), 2.7.9.3, R03595"
			     $enzyme = ucfirst $k_nums[1];
			     $enzyme =~ s/E(\d+\.\d+\.\d+\.[\d-]+)/$1/; 
			     $enzyme =~ s/(.+?)\. .*/$1/; 

			     if ( exists $deg_info{$k_nums[0]} ){
					push @{$coords_info{$enzyme}},(@{$deg_info{$k_nums[0]}},$coords);
			     }

			}else{

			    $enzyme = ucfirst $k_nums[1] unless ( $enzyme );
			    $enzyme =~ s/E(\d+\.\d+\.\d+\.[\d-]+)/$1/; 
			    $enzyme =~ s/(.+?)\. .*/$1/; 

			    # A rectangle of kegg map contains several K_numbers and these K_numbers contains both up or(and) down regulated genes,such as title="K13811 (PAPSS), K00955 (cysNC), K00956 (cysN), K00957 (cysD), K00958 (sat), 2.7.7.4, R04929" 	
			    #my ($up_count,$down_count,$normal_count) = (0,0,0);
				my ($up_count,$down_count,$normal_count)=(0,0,0);
			    foreach my $index ( 0..$#k_nums ) {
			     
			        next if $index%2 == 1;		

			        if ( exists $count_of_up_down_gene_in_K->{$k_nums[$index]} ) {
				    $count_of_up_down_gene_in_K->{$k_nums[$index]}->{up} = 0 unless $count_of_up_down_gene_in_K->{$k_nums[$index]}->{up};
				    $count_of_up_down_gene_in_K->{$k_nums[$index]}->{down} = 0 unless $count_of_up_down_gene_in_K->{$k_nums[$index]}->{down};
					$count_of_up_down_gene_in_K->{$k_nums[$index]}->{normal} = 0 unless $count_of_up_down_gene_in_K->{$k_nums[$index]}->{normal};
				    $up_count += $count_of_up_down_gene_in_K->{$k_nums[$index]}->{up} ;
				    $down_count += $count_of_up_down_gene_in_K->{$k_nums[$index]}->{down} ;
					$normal_count += $count_of_up_down_gene_in_K->{$k_nums[$index]}->{normal} ;
			 	    push @{$coords_info->{$enzyme}},$coords; 
				}
			    }

			    # judge up/down-regulation 
			    if ( $up_count > $down_count && $up_count > $normal_count) { 
				    push @{$coords_info->{$enzyme}},"up";

			    }elsif ( $up_count < $down_count && $down_count > $normal_count) { 
				    push @{$coords_info->{$enzyme}},"down";

			    }elsif($up_count==$down_count && $up_count > 0){
					push @{$coords_info->{$enzyme}},("up","down");
				}else{ 
				    push @{$coords_info->{$enzyme}},("normal");
			    }
			}	
		}
	}
	close HTML;

	# remove dupliate value
	&rmdup($coords_info);
}

sub K_number_up_down{
	my ($coords_up_down,$num,$key,$style,$div) = @_;

	#my @up_down = sort grep { $_ =~ /up|down/ } @{$coords_up_down};
	#my @all_coords = grep { $_ !~ /up|down/ } @{$coords_up_down};	
	my @up_down = sort grep { $_ =~ /up|down|normal/ } @{$coords_up_down};
	my @all_coords = grep { $_ !~ /up|down|normal/ } @{$coords_up_down};

	# A K_number only contains up or down regulated genes
	if ( @up_down == 1 ) {
	
	    foreach my $coord ( @all_coords ) {

	        $num ++	;
	        my @coords = split /,/,$coord;
                my $sig = $up_down[0]; 

                my $color = $sig =~ /up/ ? "#F00":
							$sig =~	/down/ ? "#00FF00" : "#0FF";
                $coords[0] = $coords[0] + 1;
                $coords[1] = $coords[1] + 1;
                my $width = $coords[2] - $coords[0] ;
                my $height = $coords[3] - $coords[1] ;

                $style .= "\t\t.E$num { position:absolute; left:$coords[0]px; top:$coords[1]px; background:$color; width:${width}px; height:${height}px;}\n";
                $div .= "\t\t<div class=\"E$num\" >$key</div>\n";	
	   }

	# A K_number contains both up and down regulated genes
	}else{
	 
	    foreach my $coord ( @all_coords ) {
		
		$num ++;
		my @coords = split /,/,$coord;
            	my ($sig_down,$sig_up) = @up_down;
		
	    	# first half rectangle cooords
            	my $coords_down_left  = $coords[0] + 1;
            	my $coords_down_top = $coords[1] + 1;
            	my $coords_down_width = ($coords[2] - $coords[0] - 1)/2;
            	my $coords_down_height = $coords[3] - $coords[1] - 1;
	    	$coords_down_width =~ s/\.\d+$//;
		
	    	# latter half rectangle coords
            	my $coords_up_left = $coords_down_left + $coords_down_width;
	    	my $coords_up_top = $coords_down_top;
            	my $coords_up_width = $coords[2] - $coords[0] - $coords_down_width - 1 ;
            	my $coords_up_height = $coords_down_height;

	    	# whole rectangle coords
	    	my $width = $coords[2] - $coords[0] - 1 ;
            	my $height = $coords[3] - $coords[1] - 1 ; 
		
            	$style .= "\t\t.E$num { position:absolute; left:${coords_down_left}px; top:${coords_up_top}px; width:${width}px; height:${height}px ; }\n";
		
            	$style .= "\t\t.E${num}down { position:absolute; left:${coords_down_left}px; top:${coords_up_top}px; background:#00FF00; width:${coords_down_width}px; height:${coords_down_height}px ; }\n";

            	$style .= "\t\t.E${num}up { position:absolute; left:${coords_up_left}px; top:${coords_up_top}px; background:#F00; width:${coords_up_width}px; height:${coords_up_height}px; }\n";
		
	        $div .= "\t\t<div class=\"E${num}down\" ></div>\n";	
                $div .= "\t\t<div class=\"E${num}up\" ></div>\n";	
                $div .= "\t\t<div class=\"E$num\" >$key</div>\n";	
	   }
	}

	return($num,$style,$div);
}

sub color_map_html{
	my ($mapid,$map_png,$outdir) = @_;
	
	my $style = "";
	my $div = "";
	my $num = 0;
	
	foreach my $key ( keys %coords_info ){
	
		($num,$style,$div) = &K_number_up_down(\@{$coords_info{$key}},$num,$key,$style,$div);
	}

	# html2png 
	&html2png($div,$style,$mapid,$map_png,$outdir);
}		

sub html2png{
	my ($div,$style,$mapid,$map_png,$outdir) = @_;
	
	my $out =<<"HTML";
                <html>
                <head>
                <style type="text/css">
				div{font-size:12px;  text-align:center; font-family:"Times New Roman" ;line-height:18px;letter-spacing:-1px}
				$style
				</style>
                </head>
                <body style="margin:0">
                <img src="$map_png" name="pathwayimage" name="pathwayimage" usemap="#mapdata" border="0" />
				$div
                </body>
                </html>
HTML

	open PNGHTML,">$outdir/$mapid\_png.html" or die;
	print PNGHTML "$out\n";
	close PNGHTML;

	# htm2png 
	my $cmd1 = "$wkhtmltoimage --quality 50  $outdir/$mapid\_png.html $outdir/$mapid.png";
	print "$cmd1###\n";
	
	
	my $return = system ($cmd1);
	my $count = 0;

	while ( 1 ) {	
	    $count ++;
	    if ( $return > 0 ) {
	        sleep(5);
	        $return = system($cmd1);

	    }else{
		last;
	    }

	    #max times of cycle
	    if ( $count > 3 ){
		last;
	    }
	}	

	my $cmd2 = "cp $mapdb/$mapid.html $outdir && rm $outdir/$mapid\_png.html";
	system ($cmd2);

}

sub USAGE{

	my $usage = << "USAGE";

	Program: $0

	Version: $version

	Contactor: jiangdezhi(dezhijiang\@annoroad.com)

	Usage:	perl $0 [options]

	Options:

		-h|help|?	help information
		-kegg_report	kegg enrichment file ( forced ) 
		-mapdb		directory of containing kegg map downloaded from KEGG database( forced )  
		-od		directory of containing colored map [ default: kegg_result ] 
	
	Example:
		perl $0 -kegg_report <kegg enrichment file> -mapdb <mapdb>  -od <outdir>

USAGE

	print "$usage\n";
	exit;
}
