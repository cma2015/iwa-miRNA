#!/uer/bin/perl
@ARGV=qw# ../data/all.seq.fasta-48.features #;
die "Usage:perl $0 <Input_File> \n" if(@ARGV<1);
open(INPUT,"$ARGV[0]")||die "can not open file $ARGV[0]\n";
@file1=<INPUT>;
chop(@file1);
@ARGV=qw# ../data/all.seq.fasta-104.features #;
die "Usage:perl $0 <Input_File> \n" if(@ARGV<1);
open(INPUT,"$ARGV[0]")||die "can not open file $ARGV[0]\n";
@file2=<INPUT>;
chop(@file2);

open( OUT,">../data/seq.svm.data" ) or die "result_switch.out: $!";

$flag=-1;
for($i=0;$i<@file1;$i++)
{
	@temp1=split(' ',$file1[$i]);
	@temp2=split(' ',$file2[$i]);
	$k=1;
	if($file1[$i])
	{
		print OUT ("$flag");
		for($j=0;$j<@temp1;$j++)
		{
			print OUT (" $k:$temp1[$j]");
			$k++;
		}
		for($j=0;$j<@temp2;$j++)
		{
			print OUT (" $k:$temp2[$j]");
			$k++;
		}
		print OUT ("\n");
	}
}

close OUT;
close(INPUT)||die "can not close the file:INPUT\n";

