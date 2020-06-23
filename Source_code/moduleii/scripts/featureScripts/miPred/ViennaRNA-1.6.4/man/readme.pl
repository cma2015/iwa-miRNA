#!/uer/bin/perl
@ARGV=qw# RNALfold.1 #;
die "Usage:perl $0 <Input_File> \n" if(@ARGV<1);
open(INPUT,"$ARGV[0]")||die "can not open file $ARGV[0]\n";
@file=<INPUT>;
chop(@file);

open( OUT,">RNALfold.1" ) or die "plantresult.out: $!";

  
for($j=0;$j<@file;$j++)
{
	print OUT ("$file[$j]\n");
}
close OUT;
close(INPUT)||die "can not close the file:INPUT\n";