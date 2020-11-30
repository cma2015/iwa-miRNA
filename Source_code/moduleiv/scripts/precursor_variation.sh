#!bin/sh

curdir=$1 premirna=$2; SNPs=$3; premirnaSNP=$4

## Add environment path
script_path=${curdir}/scripts

# Define the output directory
OUTDIR=${premirna%.*}/sequenceVariation

if [ ! -d ${premirna%.*} ]; then
    mkdir ${premirna%.*}
fi

if [ ! -d ${OUTDIR} ]; then
    mkdir ${OUTDIR}
fi

cd ${OUTDIR}

cp ${script_path}/sequence_complement.py sequence_complement.py
cp ${script_path}/precursor_plots.py precursor_plots.py
cp ${premirna} final_table.txt
awk 'NR>1{OFS="\t";print $1,$2-1,$2,$3,".",".",$4,$5}' ${SNPs} > SNP.bed

Rscript ${script_path}/precursor_extraction.R

bedtools intersect -wo -a pre_miRNA.bed -b SNP.bed >intersect_miSNPs.txt

python ${script_path}/precursor_SNP.py

(head -n 2 SNP.vcf && tail -n +3 SNP.vcf | sort -k 1,1 -k 2,2n ) >SNP_sort.vcf
bcftools view SNP_sort.vcf -Oz -o SNP.vcf.gz
bcftools index -f SNP.vcf.gz

if [ `grep -c "Extended_stem_loop_loc" final_table.txt` -ne '0' ];then
    awk -F"\t" 'NR>1{print ">"$2"\n"$6}' final_table.txt >pre_miRNA.fa
else
    awk -F"\t" 'NR>1{print ">"$1"\n"$3}' final_table.txt >pre_miRNA.fa
fi

python  sequence_complement.py inter_pre_miRNA.bed pre_miRNA.fa
bcftools consensus SNP.vcf.gz --fasta-ref finalPremiRNAsMod.fasta -o finalPremiRNAsMod_SNPs.fasta
python  precursor_plots.py

pdfunite $(ls -v *.pdf) ${premirnaSNP}

rm *.ps *.png *.pdf *.str
