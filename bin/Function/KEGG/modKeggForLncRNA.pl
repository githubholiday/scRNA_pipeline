my ($report) = @ARGV;

open FL,$report or die $!;
chomp(my $header = <FL>);
print "$header\tUnchange_gene\n";
while(<FL>){
	chomp;
	my (%up,%down,%de,$normal);
	my $line = $_;
	my ($de,$up,$down) = (split/\t/)[9,10,12];
	#print "$de\n";
	my @de = split/;/,$de;
	#print "@de\n";
	for (@de){
		$de{$_} = 1;
	}
	unless ($up eq "."){
		my @up = split/;/,$up;
		for (@up){
			delete $de{$_};
			$up{$_} = 1;
		}
	}
	unless ($down eq "."){
		my @down=split/;/,$down;
		for (@down){
			delete $de{$_};
			$down{$_} = 1;
		}
	}
	my @key = sort keys %de;
	$normal = join(";",@key);
	print "$line\t$normal\n";
	#print "$de\t$up\t$down\n";
	#last;
}
close FL;
