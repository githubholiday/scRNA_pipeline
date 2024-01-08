#!/usr/bin/perl -w
#use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use Cwd 'abs_path';
use File::Basename qw(basename dirname);

my ($prefix,$indir,$outdir);
GetOptions(
	"help|?" =>\&USAGE,
	"indir:s"=>\$indir,
	"prefix:s"=> \$prefix,
	"outdir:s" => \$outdir,
) or &USAGE;
#&USAGE unless ($fIn and $fOut);

my $abs_dir = abs_path($indir);
&Read_DIR($abs_dir,"DIR");

sub Read_DIR{
	my ($parent_dir,$handle) = @_;
	opendir $handle,$parent_dir or die $!; ## 打开目录
	my @files = readdir $handle;
	if ($#files > 3 or ( -d "$parent_dir/$files[2]")){ ## 目录含一个文件及以上 或只包含一个目录(因为有. 和 ..两层)
		for (@files){ ## 遍历目录
			my $file = $_;
			next if ($file =~ /^\./);
			my $file_path = $parent_dir . "/$file";
			if (-f $file_path){ ## 判断是否为文件
				&PREFIX_FILE($parent_dir,$file,$prefix); ## 如果是文件，则如此操作
				#print "file:$file\n";
			}elsif(-d $file_path){ ## 判断是否为目录
				system('mkdir -p $outdir/$file');
				&Read_DIR($file_path,$file); ## 递归读取目录
				#print "DIR:$file\n";
			}
		}
	}else{
		&PREFIX_DIR($parent_dir,$files[$#files]);
	}
	closedir $handle;
}

### 处理只含一个文件的目录
sub PREFIX_DIR{
	my ($parent_dir,$file) = @_;
	my ($new_file);
	my $dir= basename $parent_dir;
	my $dir_name = basename (dirname $parent_dir);
	$dir =~ /.*_(.*)/; ### 取目录名后缀
	my $suffix_dir = $1;
	$file =~ /(.*)\.(.*)/; ### 取文件名前缀和后缀
	my $prefix_file = $1;
	my $suffix_file = $2;
	if ($suffix_dir eq $prefix_file){ ## 如果目录后缀和文件前缀相同，则合并
		$new_file = $prefix . "\.$dir" . "\.$suffix_file";
	}else{
		$new_file = $prefix . "\.$dir" . "\.$file";
	}
	print "Link :$parent_dir/$file $outdir/$dir_name/$new_file\n";
	system("mkdir -p $outdir/$dir_name && ln -snf $parent_dir/$file $outdir/$dir_name/$new_file");
}

### 处理正常文件
sub PREFIX_FILE{
	my ($parent_dir,$file,$prefix) = @_;
	my $file_old = $parent_dir . "/$file";
	my $file_new = $prefix . "\.$file";
	#my $file_new = $parent_dir . "/$prefix" . "\.$file";
	my $dir_name = basename $parent_dir;
	print "Link :$file_old $outdir/$dir_name/$file_new\n";
	system("mkdir -p $outdir/$dir_name && ln -snf $file_old $outdir/$dir_name/$file_new");
	#rename $file_old,$file_new;
}
