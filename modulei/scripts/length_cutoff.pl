#!/usr/bin/env perl

############     length_cutoff.pl     ############
#
#	This Perl script reads sequence files in FASTA format
#	and sorts the sequences into three output files according
#	to the applied length cutoff:
#	- sequence length ok (sequences_ok.fas)
#	- sequence too short (sequences_too_short.fas)
#	- sequence too long (sequences_too_long.fas)
#
#	The script is part of the "NGS tools for the novice"
#	authored by David Rosenkranz, Institue of Anthropology
#	Johannes Gutenberg University Mainz, Germany
#
#	Contact: rosenkrd@uni-mainz.de



############     HOW TO USE     ############
#
#	You can pass file names of the files to be processed and
#	the values of the length cufoff as arguments to the script.
#
#	Input files have to be passed to the scripts via -i:
#	-i file1.fas -i file2.fas
#
#	Cutoff values have to be passed to the script via -min and
#	-max respectively:
#	-min 20 -max 70
#
#	A bare digit (0 or 1) will tell the script how to treat
#	sequences that exceed the maximum length. By default
#	it is set to 0.
#	-0 = remove sequences (write to output file sequences_too_long.fas)
#	-1 = cut the end of the sequence and save to file sequences_ok.fas
#
#	For example you can type the following command:
#	perl length_cutoff.pl -i file1.fas -min 30 -max 60 -1
#
#	If you do not want to enter each file name seperately, you
#	can provide a file that contains a list all file names (one
#	file name per line). Pass the file name to the script via -I:
#	-I list_of_files.txt
#
#	Multiple files and combinations of all kinds of arguments are
#	allowed:
#	perl length_cutoff.pl -i file1.fas -I list_of_files.txt -min 30 -max 60 -1




@input_files=();
$min_length=0;
$max_length=999999999999;
$treat_long=0;
$|=1;

###   CHECK COMMAND LINE ARGUMENTS   ###
if(@ARGV==0)
	{
	print"No arguments passed to the script!\nIf you entered arguments try the following command:\nperl sort_by_tag.pl -argument1 #argument2 ...\n\n";
	exit;
	}

$argv="";
foreach(@ARGV)
	{
	$argv.=$_;
	}
@arguments=split('-',$argv);

foreach(@arguments)
	{
	if($_=~/^ *i/)
		{
		$_=~s/^ *i//;
		$_=~s/ //g;
		push(@input_files,$_);
		}
	elsif($_=~/^ *I/)
		{
		$_=~s/^ *I//;
		$_=~s/ //g;
		open(FILES_IN,"$_");
		while(<FILES_IN>)
			{
			unless($_=~/^\s*$/)
				{
				$_=~s/\s//sg;
				push(@input_files,$_);
				}
			}
		}
	elsif($_=~/^ *min/)
		{
		$_=~s/^ *min//;
		$_=~s/ //g;
		$min_length=$_;
		unless($min_length=~/^\d+$/)
			{
			print"Minimum sequence length has to be numerical!\n";
			exit;
			}
		}
	elsif($_=~/^ *max/)
		{
		$_=~s/^ *max//;
		$_=~s/ //g;
		$max_length=$_;
		unless($max_length=~/^\d+$/)
			{
			print"Maximum sequence length has to be numerical!\n";
			exit;
			}
		}
	elsif($_=~/^ *name/)
		{
		$_=~s/^ *name//;
		$_=~s/ //g;
		$nameIndex=$_;
		}
	elsif($_=~/^ *[01]/)
		{
		$_=~s/ *//g;
		$_=~s/ //g;
		$treat_long=$_;
		}
	elsif($_!~/^\s*$/)
		{
		print"Don't know how to treat argument $_!\nIt will be ignored.\n\n";
		}
	}
unless($min_length<=$max_length)
	{
	print"Maximum sequence length must be >= minimum seuquence length!\n";
	exit;
	}
if(@input_files==0)
	{
	print"No input file specified!\n";
	exit;
	}

###   PRINT ARGUMENTS   ###
print"The following files will be processed:\n";
foreach(@input_files)
	{
	if(-e $_)
		{
		print"$_\n";
		push(@input_files_ok,$_);
		}
	else
		{
		print"could not find file: $_. It will be ignored.\n";
		}
	}
print"Minimum seuquence length: $min_length\nMaximum seuquence length: $max_length\n";
if($treat_long==1)
	{
	print"Sequences exceeding maximum length will be clipped to proper size.\n\n";
	}

###   START   ###
open(OUT_OK,">$nameIndex");
unless($treat_long==1)
	{
	open(OUT_LONG,">$nameIndex_long.fa");
	}
open(OUT_SHORT,">$nameIndex_short.fa");
$seqs_ok=0;
$seqs_long=0;
$seqs_short=0;
foreach$file(@input_files_ok)
	{
	print"processing $file";
	open(IN,$file);
	while(<IN>)
		{
		if($_=~/^>/)
			{
			$title=$_;
			}
		elsif($_!~/^\s*$/)
			{
			chomp$_;
			if(length$_>=$min_length)
				{
				if(length$_<=$max_length)
					{
					$seqs_ok++;
					print OUT_OK "$title$_\n";
					}
				elsif($treat_long==0)
					{
					$seqs_long++;
					print OUT_LONG "$title$_\n";
					}
				elsif($treat_long==1)
					{
					$_=substr($_,0,$max_length);
					$seqs_long++;
					print OUT_OK "$title$_\n";
					}
				}
			else
				{
				$seqs_short++;
				print OUT_SHORT "$title$_\n";
				}
			}
		}
	close IN;
	print" done.\n";
	}
close OUT_OK;
close OUT_LONG;
close OUT_SHORT;
print"Sequences with proper size:\t$seqs_ok\nSequences too long:\t\t$seqs_long\nSequences too short:\t\t$seqs_short\n\n";
exit;
