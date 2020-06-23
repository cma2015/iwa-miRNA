#!/uer/bin/perl
@ARGV=qw# ../data/seq.fasta.fold #;
die "Usage:perl $0 <Input_File> \n" if(@ARGV<1);
open(INPUT,"$ARGV[0]")||die "can not open file $ARGV[0]\n";
@file=<INPUT>;
chop(@file);
open( OUT,">../data/seq_mat.fasta" ) or die "result_switch.out: $!";

for($i=0;$i<@file;$i+=3)
{
	@temp1=split('',$file[$i+1]);
	@temp2=split('',$file[$i+2]);
	$leng=@temp1;
	for ($j1=0,$j2=$leng;$j2-$j1>60;)
	{
		print OUT ("$file[$i]\n");
		for($k=$j1;$k<$j2;$k++)
		{
			print OUT ("$temp1[$k]");
		}
		print OUT ("\n");
		$count=0;
		if($temp2[$j1] eq '(')
		{
			$count++;
		}
		if($temp2[$j1+1] eq '(')
		{
			$count++;
		}
		if($temp2[$j1+2] eq '(')
		{
			$count++;
		}
		$j1+=3;
		while($count>0)
		{
			if($temp2[$j2-1] eq ')')
			{
				$count--;
			}
			$j2--;
		}
	}
}
close OUT;
close(INPUT)||die "can not close the file:INPUT\n";
