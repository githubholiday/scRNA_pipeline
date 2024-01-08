#!/usr/bin/perl
use Encode;
use diagnostics;
use Mail::Sender;
use Getopt::Long;
#use utf8;
my (%opts);
GetOptions(\%opts,"input:s","flag:s","subject:s" ,"usrmail:s" ,"help");

sub help{
        print STDERR <<USEAGE;
        perl $0 -input <txt> -flag <auto> -subject <Email Subject> -usrmail <Mail Receiver> 
                -input                        邮件正文(文本文件或者数据库) example: "-input QUALITY,VAULE"
                -flag                         <yes/no/[ABCD]>  default:yes（yes:发信；no:不发；[ABCD]:判断发邮件）
                -subject                      邮件主题
                -usrmail                      邮件配置文件

USEAGE
}
unless($opts{input} and $opts{subject} and $opts{usrmail} ){
        help;
        exit;
}

my $auto=$opts{flag};
$auto||='A';
my $flag||='A';

my ($file,$value)=(split(/,/,$opts{input}))[0,1];

my $subject=$opts{subject};
#my $usrmail=$opts{usrmail};
$subject=encode("gbk", decode("utf-8", $subject));

#print "$file\t$subject\n";

my $type=`file $file`;
my $x;

if( $type =~ 'text'){
    $x=&Txt($file,$value);
}elsif($type =~ 'SQLite' ){
        $x=&Sqlite($file);
}else{
    die "文件不对!";
}

#print "STANDARD:$auto\nRESULT:$flag\n";

if( $auto =~ 'no'){
    #print 'Good Bye!';
    exit;
}

if( $auto =~ 'yes'){
    &Mail;
}


if( $auto le $flag ){
    &Mail;
}

##############发邮邮件###################

sub Mail{

my($Addressor,$Password,$Receiver,$Copy);
#open IN,"$ENV{HOME}/.email/.email.txt" or die "PLEASE DEFINED YOU EMAIL CONFIG!!!";
open IN,$opts{usrmail} or die "PLEASE DEFINED YOU EMAIL CONFIG!!!";
while(<IN>){
	chomp;
	if(/Addressor/){
	$Addressor=(split(/=/,$_,2))[1];
	$Addressor=~s/\s//g;
	}
	if(/Password/){
	$Password=(split(/=/,$_,2))[1];
	$Password=~s/\s//g;
	}
	if(/Receiver/){
	$Receiver=(split(/=/,$_,2))[1]; 
	$Receiver=~s/;/,/g;
	$Receiver=~s/\s//g;
	}
	if(/Copy/){
	$Copy=(split(/=/,$_,2))[1];
	$Copy=~s/;/,/g;   
	$Copy=~s/\s//g;
	}
}
close IN;

my $filegbk=encode("gbk", decode("utf-8", $file));
$filegbk="$filegbk,$value";
eval {
(new Mail::Sender)
#->MailMsg({
->OpenMultipart({
    smtp => 'smtp.exmail.qq.com',#hostname or ip 
    from => $Addressor, 
#    fake_from => 'tecent@qq.com',
    priority => '1',
    confirm => 'delivery, reading',
#    debug => 'mailbug.txt', 
#    to=>$usrmail,
	to => $Receiver,
#    cc => $Copy,
    subject => $subject,
    auth => 'LOGIN',#验证方式
#    auth_encoded => 'true',
#    authid => 'c2VxeXVhbkBhbm5vcm9hZC5jb20=',
#    authid => 'seqyuan@annoroad.com',
    authid => $Addressor,
    authpwd => $Password,
#    authpwd => 'U2VxeXVhbjk5OQ==',
})

    ->Part({ctype => 'text/html', disposition => 'NONE',msg => $x})


    ->Attach({
    disposition => "inline;filename=$filegbk;\r\nContent-ID: <img1>",
    file => $filegbk
    })
    ->Close();
} or print "Error sending mail: $Mail::Sender::Error\n";

#print"Mail sent OK.\n";

}


sub Txt{
    my $i=shift;
    my $ii=shift;
    my $s;
    my $re='';
    open IN,$i or die;
    while(<IN>){
        chomp;
        if(/\b[BCD]\b/){
		$re.="<p>$_ Failed</p>\n";
		}
        my @tab=split(/\t/,$_);
        my $autochar=(split(/\t/,$_,2))[1];
        my @tt=split(/\t/,$autochar);
        $s.='<tr>';
        for my $h (@tab){
            for my $k(@tt){
                if($k=~ /\b([ABCD])\b/ && $flag lt $k){
                    $flag=$k;
                }
            }
				if($h eq 'A'){
				$h='通过';
				}elsif($h eq 'B'){
				$h='让步';
				}elsif($h eq 'C'){
				$h='不合格';
				}elsif($h eq 'D'){
				$h='终止';
				}
                $s.=qq(<td align="center" >$h</td>\n);
            }
        $s.='</tr>';
    }

    close IN;
	my $ss;
	open IN,$ii or die;
	while(<IN>){
	chomp;
	my @tab=split(/\t/,$_);
	$ss.='<tr>';
	for my $h (@tab){
		$ss.=qq(<td align="center" >$h</td>\n);
	}
	$ss.='</tr>';
	}
	close IN;
my $txt='
<html >
 <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
 </head>
<body>
<table border="1" align="left" cellspacing="0">
 <caption><b>QUALITY</b></caption> 
     '.$s.
' 
<!--/table-->
<b>'
.$re.
'
<b>
</table>
<br />
<table border="1" align="left" cellspacing="0">
 <caption><b>VALUE</b></caption>'
.$ss.
'</table> 
</body>
</html>
';
    return $txt;
}

sub Sqlite{
    my $database=shift;
    my $s;
    my $c=`sqlite3 $database '.table' `;
    my @get=split(/\s+/,$c);
    my @t;
    for my $get (@get){
        if($get=~/QUALITY/){
            push @t,$get;
        }
    }
    for my $t (@t){
    my $sn=`sqlite3  -header  -csv $database 'select * from $t' `;
    $s.='<div><p>&nbsp</p>';
    $s.='<table border="1" align="left" cellspacing="0">';
    $s.='<caption ><b>'.$t.'</b></caption>';
    my @sn=split(/\n/,$sn);

    for my $i(@sn){
    my @tab=split(/,/,$i);
    my $autochar=(split(/,/,$i,2))[1];
    my @tt=split(/,/,$autochar);
    for my $k(@tt){
       if($k=~ /\b([ABCD])\b/ && $flag lt $k){
           $flag=$k;
       }
    }
    $s.='<tr>';
    for my $h (@tab){
            $s.=qq(<td align="center" >$h</td>\n);
    }
    $s.='</tr>';
    }
    $s.='</table>';
    $s.='<p>&nbsp</p><p>&nbsp</p></div>';
    }

    my $sqtxt='
<html >
 <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
 </head>
<body>
     '.$s.
' 
</body>
</html>
';

    return $sqtxt;
}

