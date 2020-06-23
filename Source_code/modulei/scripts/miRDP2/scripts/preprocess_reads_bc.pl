#!/bin/perl

use strict;
use Getopt::Std;


######################################### USAGE ################################

my $usage =
"$0 reads_file rfam_aln miR_aln threshold processed_reads filtered_reads

This script is used to preprocess and filter original reads file of miRDeep-P.
The file, reads_file, is the fasta format unique reads file, with id formatted as '>reads00001_x123'
The file, rfam_aln, is the default output of bowtie-mapping of all reads against ncRNA seq (tRNA+rRNA+snRNA) from Rfam
The file, miR_aln, is the default output of bowtie-mapping of all reads against known plant miRNA mature seq from miRbase
The number, threshold, is the rpm threshold for reads to be retrieved for precursor excision
The file, processed_reads, is the output path of processed reads file for signature preparation
The file, filtered_reads, is the output path of filtered reads file for precursor excision
";


####################################### INPUT FILES ############################

my $reads_file=shift or die $usage;
my $rfam_aln=shift or die $usage;
my $miR_aln=shift or die $usage;
my $threshold =shift or die $usage;

my $processed_file=shift or die $usage;
my $filtered_file =shift or die $usage;
my $total_reads_file =shift or die $usage;


##################################### GLOBAL VARIBALES #########################

my $reads_len_min=19;
my $reads_len_max=24;
my $total_reads=0;

my %reads_hash=();
my %ncRNA_hash=();
my %miR_hash=();


######################################### MAIN ################################# 


parse_aln($rfam_aln,\%ncRNA_hash);

parse_aln($miR_aln,\%miR_hash);

parse_reads($reads_file,\%reads_hash);


open(PROCESSED,">",$processed_file);
open(FILTERED,">",$filtered_file);
open(TOTALREADS,">",$total_reads_file);

# filter reads by sequence similarity to miRNA/ncRNA & reads rpm
foreach my $id(sort keys %reads_hash)
{
	#print all non-ncRNA reads for signature preparaion
	print PROCESSED ">$id\n$reads_hash{$id}\n";
	
	#retrieve reads match known miRNA mature seq and exceeded rpm threshold
	$id=~/_x(\d+)/;
	if($miR_hash{$id} eq "T" || $1/$total_reads*1000000 >= $threshold){print FILTERED ">$id\n$reads_hash{$id}\n";}
}
print TOTALREADS "$total_reads";



######################################### SUBROUTINES ################################# 

# parse aln file, extract reads ids correlated with rfam ncRNA or known miRNA
sub parse_aln
{
	my($file,$in_hash)=@_;
	my $filehandle="ALN_HANDLE";
	open(ALNHANDLE,"$file") or die "open file failed: $!\n";
	while(<ALNHANDLE>){my @a=split "\t"; $$in_hash{$a[0]}="T";}
	close(ALNHANDLE);
}

# filter reads & count total reads number
sub parse_reads
{
	my($file,$in_hash)=@_; 
	my $id=""; my $seq="";
	
	open(READSHANDLE,"$file") or die "open file failed: $!\n";
	while(<READSHANDLE>)
	{chomp;s/\r//;s/^>//;$id=$_;$seq=<READSHANDLE>;chomp($seq);$seq=~s/\r$//; 
		$id=~/_x(\d+)/;$total_reads+=$1;
	
	 #discard all ncRNA reads
	 if($ncRNA_hash{$id} eq "T"){next;}
	 if(length($seq)>=$reads_len_min && length($seq)<=$reads_len_max)
	 {$$in_hash{$id}=$seq; 	#$id=~/_x(\d+)/; $total_reads+=$1;
	 } 
	}
	close(READSHANDLE);
}

