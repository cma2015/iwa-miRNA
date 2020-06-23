#!/uer/bin/perl
#@ARGV=qw# ../data/seq.fasta.fold #;
die "Usage:perl $0 <Input_File> \n" if(@ARGV<1);
open(INPUT,"$ARGV[0]")||die "can not open file $ARGV[0]\n";
@file1=<INPUT>;
chop(@file1);

open( OUT,">../data/seq.fasta.fold1" ) or die "result_switch.out: $!";

for($i=0;$i<@file1;$i+=3)
{
	@temp=split(' ',$file1[$i+2]);
	@temp1=split('\(',$temp[1]);
	$leng=@temp;
	print ("$leng\n");
	if($leng == 3)
	{
		@temp2=split('\)',$temp[2]);
	}
	else
	{
		@temp2=split('\)',$temp1[1]);
	}
	$leng=length($file1[$i+1]);
	print OUT ("$file1[$i]\n$file1[$i+1]\n$temp[0]\n$temp2[0]\n$leng\n");
}

close OUT;
close(INPUT)||die "can not close the file:INPUT\n";

