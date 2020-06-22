#!/bin/bash

# The main pipeline of miRDP2


usage="$0 <OPTIONS>

PLEASE REFER TO MIRDP2 MANUAL ALSO.

  OPTIONS:
    NECESSARY:
    -g/--genome <file>        reference genome file in fasta format
    -i/--index <path>/prefix  bowtie index file corresponding to reference genome
    -f/--fasta <file>         formatted sRNA-seq file in fasta format
    -o/--output <path>        path to output folder

    OPTIONAL:
    -L/--locate <int>         THRESHOLD of different locations one read can map to, default is 15
    -M/--mismatch <int>       THRESHOLD of allowed mismatches for bowtie mapping, default is 0
    -R/--rpm <int>            THRESHOLD of reads rpm when filtering, default is 10
    -p/--thread <int>		  number of thread used for RNAfold, default is 1
	--large-index             use this option when using large bowtie index, e.g. when using *.ebtwl files
    -h/--help                 print this help
"


#####################  default parameters  #####################
locate=15
mismatch=0
rpm=10
large=""
thread=1


#####################  get input  #####################
while true
do
  if [[ "x$1" == "x" ]]; then break;
  fi;

  case "$1" in
    -d|--dir)  dir=$2; shift 2;;	#input directory
	-g|--genome)  genome=$2; shift 2;;	#input genome file
    -i|--index)   bowtie_index=$2; shift 2;;	#input index file
    -f|--fasta)   seq=$2; shift 2;;	#input fastq seq file
    -o|--output)  results_folder=$2; shift 2;;	#output folder; default is `.'
    -L|--locate)  locate=$2; shift 2;;	#multi_location threshold
    -M|--mismatch) mismatch=$2; shift 2;;	#mismatch
    -R|--rpm)      rpm=$2; shift 2;;	#RPM threshold
    -p|--thread)      thread=$2; shift 2;;	#thread number
	--large-index)    large="--large-index"; shift;;
    -h|--help)    echo -e "$usage"; exit 1;;
    *)            echo -e "$usage"; exit 1;;
  esac
done



#####################  change variables  #####################
filename=`echo $seq | perl -e '$a=<>;@b=split "/",$a;$b[-1]=~/(\S+)\.\w+$/;print $1'`	#seq file prefix

src=$dir #${0%%/miRDP2-v1.1.1_pipeline.bash}

mirdp=$src/scripts			#miRDP script folder
len=$locate				#allowed maximum mapping site number
var=$mismatch				#allowed mismatch number
threshold=$rpm		#threshold for filter reads (rpm)

filename=$filename-$len-$var-$threshold

#####################  begin pipeline  #####################
mkdir -p $results_folder/$filename
echo "begin..." >> $results_folder/${filename}/progress_log
date +%Y-%m-%d=%X >> $results_folder/${filename}/progress_log


#filter reads -- mapping to ncRNA seq (rRNA, tRNA, etc)
#bowtie -v 0 $mirdp/index/rfam_index -f $input_folder/${1}.fa > $results_folder/${filename}/rfam_reads.aln
bowtie -v 0 $mirdp/index/rfam_index -f $seq > $results_folder/${filename}/rfam_reads.aln 2>> $results_folder/${filename}/script_err

#filter reads -- mapping to known miR mature seq
#bowtie-build -f $mirdp/plant_isoform.fa $mirdp/index/mature_index
bowtie -v 1 $mirdp/index/mature_index -f $seq > $results_folder/${filename}/known_miR.aln 2>> $results_folder/${filename}/script_err

#filter reads -- extract reads with high count number
perl $mirdp/preprocess_reads.pl $seq $results_folder/${filename}/rfam_reads.aln $results_folder/${filename}/known_miR.aln $threshold $results_folder/${filename}/${filename}.fa $results_folder/${filename}/${filename}-processed.fa $results_folder/${filename}/${filename}.total_reads 2>> $results_folder/${filename}/${filename}_err


#mapping filtered reads
bowtie -a -v $var $large $bowtie_index -f $results_folder/${filename}/${filename}-processed.fa > $results_folder/${filename}/${filename}_processed.aln 2>> $results_folder/${filename}/script_err


#filter reads -- ignore reads mapping to too many sites
perl $mirdp/convert_bowtie_to_blast.pl $results_folder/${filename}/${filename}_processed.aln $results_folder/${filename}/${filename}.fa $genome > $results_folder/${filename}/${filename}-processed.bst 2>>$results_folder/${filename}/${filename}_err
perl $mirdp/filter_alignments.pl $results_folder/${filename}/${filename}-processed.bst -c $len > $results_folder/${filename}/${filename}-processed_filter${len}.bst 2>>$results_folder/${filename}/${filename}_err

echo "finish filtering" >> $results_folder/${filename}/progress_log

#excise precursor candidates and predict secondary strucure
perl $mirdp/excise_candidate.pl $genome $results_folder/${filename}/${filename}-processed_filter${len}.bst 250 > $results_folder/${filename}/${filename}_precursors.fa 2>>$results_folder/${filename}/${filename}_err;

if [ "$thread" -gt 1 ]
	then
		RNAfold --noPS -j$thread < $results_folder/${filename}/${filename}_precursors.fa > $results_folder/${filename}/${filename}_structures 2>>$results_folder/${filename}/${filename}_err;
	else
		RNAfold --noPS < $results_folder/${filename}/${filename}_precursors.fa > $results_folder/${filename}/${filename}_structures 2>>$results_folder/${filename}/${filename}_err;
fi

echo "finish folding candidate precursors" >> $results_folder/${filename}/progress_log


#extract reads with no ncRNA for signature preparation
#perl -e '$reads=shift;$idf=shift;open(READS,$reads);open(ID,$idf); while(<ID>){chomp;$hash{">".$_}="T";} while(<READS>){chomp;$id=$_;$seq=<READS>;chomp($seq); if($hash{$id} eq "T"){next} else{print "$id\n$seq\n";} }' $results_folder/${filename}/${filename}.fa $results_folder/${filename}/rfam-id_file > $results_folder/${filename}/${filename}-rm_nc.fa
bowtie -a -v $var $large $bowtie_index -f $results_folder/${filename}/${filename}.fa > $results_folder/${filename}/${filename}.aln 2>> $results_folder/${filename}/script_err
perl $mirdp/convert_bowtie_to_blast.pl $results_folder/${filename}/${filename}.aln $results_folder/${filename}/${filename}.fa $genome > $results_folder/${filename}/${filename}.bst 2>>$results_folder/${filename}/${filename}_err
perl $mirdp/filter_alignments.pl $results_folder/${filename}/${filename}.bst -c $len > $results_folder/${filename}/${filename}_filter${len}.bst 2>>$results_folder/${filename}/${filename}_err
perl $mirdp/filter_alignments.pl $results_folder/${filename}/${filename}_filter${len}.bst -b $results_folder/${filename}/${filename}.fa > $results_folder/${filename}/${filename}_filtered.fa 2>>$results_folder/${filename}/${filename}_err


#prepare reads signature file
mkdir $results_folder/${filename}/index
bowtie-build -f $results_folder/${filename}/${filename}_precursors.fa $results_folder/${filename}/index/${filename}_precursors >> $results_folder/${filename}/script_log 2>> $results_folder/${filename}/script_err
bowtie -a -v $var $results_folder/${filename}/index/${filename}_precursors -f $results_folder/${filename}/${filename}_filtered.fa > $results_folder/${filename}/${filename}_precursors.aln 2>> $results_folder/${filename}/script_err
perl $mirdp/convert_bowtie_to_blast.pl $results_folder/${filename}/${filename}_precursors.aln $results_folder/${filename}/${filename}_filtered.fa $results_folder/${filename}/${filename}_precursors.fa > $results_folder/${filename}/${filename}_precursors.bst 2>>$results_folder/${filename}/${filename}_err
sort +3 -25 $results_folder/${filename}/${filename}_precursors.bst > $results_folder/${filename}/${filename}_signatures 2>>$results_folder/${filename}/${filename}_err
echo "finish assigning candidate signatures" >> $results_folder/${filename}/progress_log


#miRDP core algorithm
perl $mirdp/mod-miRDP.pl $results_folder/${filename}/${filename}_signatures $results_folder/${filename}/${filename}_structures > $results_folder/${filename}/${filename}_predictions 2>>$results_folder/${filename}/${filename}_err
echo "finish core algorithm" >> $results_folder/${filename}/progress_log



#calculate chr_length
perl -e '%hash=();@out=();   $name=<>; chomp($name);s/\s+$//; $seq=();  while(<>)  {chomp;  s/\s+$//;  if(m/^[A-Za-z\*]+$/){$seq.=$_;}   elsif(m/^>/){$len=length($seq);$name=~/^>(\S+)\s*/;print "$1\t$len\n"; $name=$_; $seq="";} elsif(m/^$/){next;}  else{print "err! $_\n";} } $len=length($seq);$name=~/^>(\S+)\s*/;print "$1\t$len\n";' $genome > $results_folder/${filename}/chr_length

#remove redundant and filter by plant-specific criteria
perl $mirdp/mod-rm_redundant_meet_plant.pl $results_folder/${filename}/chr_length $results_folder/${filename}/${filename}_precursors.fa $results_folder/${filename}/${filename}_predictions $results_folder/${filename}/${filename}.total_reads $results_folder/${filename}/${filename}_nr_prediction $results_folder/${filename}/${filename}_filter_P_prediction 2>>$results_folder/${filename}/${filename}_err
echo "finish filtering plant criteria" >> $results_folder/${filename}/progress_log


#mkdir $results_folder/${filename}/tmp_files
#mv $results_folder/${filename}/* $results_folder/${filename}/tmp_files 2> /dev/null
#cp $results_folder/${filename}/tmp_files/${filename}_filter_P_prediction $results_folder/${filename}

perl -e ' while(<>){chomp;@tmp=split "\t",$_; $tmp[4]=~m/^(\d+)\.\.(\d+)$/;$mature_beg=$1; $tmp[5]=~m/^(\d+)\.\.(\d+)$/;$pre_beg=$1; if($tmp[1] eq "+"){$pre_beg-=2;$mature_beg-=2;} elsif($tmp[1] eq "-"){$pre_beg+=1;$mature_beg+=1;} $mature_len=length($tmp[6]);$pre_len=length($tmp[7]); $mature_end=$mature_bed+$mature_len;$pre_end=$pre_beg+$pre_len; $sign=""; if(($mature_beg==$pre_beg && $tmp[1] eq "+")||($mature_end==$pre_end && $tmp[1] eq "-") ){$sign="5";}elsif(($mature_beg==$pre_beg && $tmp[1] eq "-")||($mature_end==$pre_end && $tmp[1] eq "+")){$sign="3";}else{$sign="A";}  print "$tmp[0]\t$pre_beg\t$pre_end\t$tmp[3]\t$sign\t$tmp[1]\n";} ' $results_folder/${filename}/${filename}_filter_P_prediction > $results_folder/${filename}/${filename}_filter_P_prediction.bed 2>>$results_folder/${filename}/${filename}_err



echo "end..." >> $results_folder/${filename}/progress_log
date +%Y-%m-%d=%X >> $results_folder/${filename}/progress_log
