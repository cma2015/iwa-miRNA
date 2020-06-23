#!bin/sh

codePath=$1 premirnaBed=$2; premirnaSeq=$3; matureLoc=$4; SNPs=$5 premirnaSNP=$6

## Create folder
parentPath=`dirname ${codePath}`
fullname=${2##*/}
fileName=${fullname%.*}
mkdir -p ${parentPath}/tmp/${fileName}/SNPs

## copy file to file folder
cd ${parentPath}/tmp/${fileName}/SNPs
cp ${codePath}/sequenceComplement.py ./
cp ${codePath}/Hairpin_plots.py ./
cp ${matureLoc} mature_location.txt
cp ${premirnaBed} premirLoc.bed

zcat ${SNPs} | grep -P "##|TSA=SNV" | bedtools intersect -wao  -a ${premirnaBed} -b - | \
awk 'BEGIN{OFS="\t";print "##fileformat=VCFv4.2""\n""#CHROM","POS","ID","REF","ALT"}{gsub(/T/,"U",$10);gsub(/T/,"U",$11);if($6=="+"){print $4,$8-$2,$9,$10,
$11}else{print $4,$3-$8+1,$9,$10,$11}}' - > SNP.vcf
(head -n 2 SNP.vcf && tail -n +3 SNP.vcf | sort -k 1,1 -k 2,2n ) >SNP_sort.vcf
bcftools view SNP_sort.vcf -Oz -o SNP.vcf.gz
bcftools index -f SNP.vcf.gz
#bcftools consensus SNP.vcf.gz --fasta-ref finalPremiRNAs.fasta -o finalPremiRNAs_SNPs.fasta

python sequenceComplement.py ${premirnaBed} ${premirnaSeq}
bcftools consensus SNP.vcf.gz --fasta-ref finalPremiRNAsMod.fasta -o finalPremiRNAsMod_SNPs.fasta
python Hairpin_plots.py 

pdfunite $(ls -v *.pdf) ${premirnaSNP}
