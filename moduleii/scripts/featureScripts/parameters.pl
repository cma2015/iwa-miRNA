#!/uer/bin/perl
#@ARGV=qw#  ../data/seq.fasta.fold1 #;
die "Usage:perl $0 <Input_File> \n" if(@ARGV<1);
open(INPUT,"$ARGV[0]")||die "can not open file $ARGV[0]\n";
@file=<INPUT>;
chop(@file);

open( OUT,">../data/all.seq.fasta-104.features" ) or die "result_switch.out: $!";

for($i=0;$i<@file;$i+=5)
{
	$leng=$file[$i+4];
	$mfe=$file[$i+3];
	@temp1=split('',$file[$i+1]);
	@temp2=split('',$file[$i+2]);

	$gc_stem=0;
	$tot_base=0;
	for($j=0;$j<$leng;$j++)
	{
		if($temp1[$j] eq 'C' || $temp1[$j] eq 'G')
		{
			$gc_stem++;
		}
		if($temp2[$j] eq '(')
		{
			$tot_base++;
		}
	}
	$mfe5=$mfe/($gc_stem/$leng);
	$mfe6=$mfe/$tot_base;

	$mis_num=0;
	for($j=0;$j<$leng-21;$j++)
	{
		for($k=0;$k<21;$k++)
		{
			if($temp2[$j+k] eq '.')	
			{
				$mis_num++;
			}
		}
	}
	$avg_mis_num=$mis_num/($leng-20);

	$gc_begin=0;
	$gc_end=0;
	$mis_num_begin=0;
	$mis_num_end=0;
	for($j=0;$j<21;$j++)
	{
		if($temp1[$j] eq 'C' || $temp1[$j] eq 'G')
		{
			$gc_begin++;
		}
		if($temp2[$j] eq '.')
		{
			$mis_num_begin++;
		}
	}
	for($j=$leng-21;$j<$leng;$j++)
	{
		if($temp1[$j] eq 'C' || $temp1[$j] eq 'G')
		{
			$gc_end++;
		}
		if($temp2[$j] eq '.')
		{
			$mis_num_end++;
		}
	}
	if($gc_begin==0)
	{
		$gc_begin=0.1;
	}
	$mfe7=$mfe/($gc_begin/21);
	if($gc_end==0)
	{
		$gc_end=0.1;
	}
	$mfe8=$mfe/($gc_end/21);
	$mfe9=$mfe/$avg_mis_num;
	
	@array_begin;
	@array_end;
	@array;
	for($j=0;$j<4;$j++)
	{
		for($k=0;$k<8;$k++)
		{
			$array_begin[$j][$k]=0;
			$array_end[$j][$k]=0;
			$array[$j][$k]=0;
		}
	}
	for($j=0;$j<21;$j++)
	{
		$temp_num1=0;
		$temp_num2=0;
		if($temp1[$j] eq 'C')
		{
			$temp_num1+=1;
		}
		if($temp1[$j] eq 'G')
		{
			$temp_num1+=2;
		}
		if($temp1[$j] eq 'U')
		{
			$temp_num1+=3;
		}
		if($temp2[$j] eq '(' || $temp2[$j] eq ')')
		{
			$temp_num2+=2;
		}
		if($temp2[$j+1] eq '(' || $temp2[$j+1] eq ')')
		{
			$temp_num2+=4;
		}
		if($j!=0 && ($temp2[$j-1] eq '(' || $temp2[$j-1] eq ')'))
		{
			$temp_num2++;
		}
		$array_begin[$temp_num1][$temp_num2]++;
	}
	for($j=$leng-21;$j<$leng;$j++)
	{
		$temp_num1=0;
		$temp_num2=0;
		if($temp1[$j] eq 'C')
		{
			$temp_num1+=1;
		}
		if($temp1[$j] eq 'G')
		{
			$temp_num1+=2;
		}
		if($temp1[$j] eq 'U')
		{
			$temp_num1+=3;
		}
		if($temp2[$j] eq '(' || $temp2[$j] eq ')')
		{
			$temp_num2+=2;
		}
		if($temp2[$j-1] eq '(' || $temp2[$j-1] eq ')')
		{
			$temp_num2+=4;
		}
		if(($j+1!=$leng) &&( $temp2[$j+1] eq '(' || $temp2[$j+1] eq ')'))
		{
			$temp_num2++;
		}
		$array_end[$temp_num1][$temp_num2]++;
	}
	for($j=0;$j<$leng;$j++)
	{
		$temp_num1=0;
		$temp_num2=0;
		if($temp1[$j] eq 'C')
		{
			$temp_num1+=1;
		}
		if($temp1[$j] eq 'G')
		{
			$temp_num1+=2;
		}
		if($temp1[$j] eq 'U')
		{
			$temp_num1+=3;
		}
		if($temp2[$j] eq '(' || $temp2[$j] eq ')')
		{
			$temp_num2+=2;
		}
		if($j-1!=0 && ($temp2[$j-1] eq '(' || $temp2[$j-1] eq ')'))
		{
			$temp_num2+=4;
		}
		if($j+1!=$leng && ($temp2[$j+1] eq '(' || $temp2[$j+1] eq ')'))
		{
			$temp_num2++;
		}
		$array[$temp_num1][$temp_num2]++;
	}
	print OUT ("$mfe5 $mfe6 $avg_mis_num ");
	print OUT ("$mfe7 $mfe8 $mfe9 $mis_num_begin $mis_num_end");
	for($j=0;$j<4;$j++)
	{
		for($k=0;$k<8;$k++)
		{
			print OUT (" $array_begin[$j][$k]");
		}
	}
	for($j=0;$j<4;$j++)
	{
		for($k=0;$k<8;$k++)
		{
			print OUT (" $array_end[$j][$k]");
		}
	}
	for($j=0;$j<4;$j++)
	{
		for($k=0;$k<8;$k++)
		{
			$array[$j][$k]=$array[$j][$k]/$leng;
			print OUT (" $array[$j][$k]");
		}
	}
	print OUT ("\n");
}

close OUT;
close(INPUT)||die "can not close the file:INPUT\n";


