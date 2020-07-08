#!bin/sh

curdir=$1 premirna=$2; SNPs=$3; premirnaSNP=$4

## Add environment path
script_path=${curdir}/scripts

# Define the output directory
out_dir=${curdir}/../tmp/sequenceVariation_`date +%s%N | cut -c6-13`

if [ ! -d ${out_dir} ]; then
    mkdir ${out_dir}
fi

cd ${out_dir}

cp ${script_path}/sequenceComplement.py sequenceComplement.py
cp ${script_path}/Hairpin_plots.py Hairpin_plots.py
cp ${premirna} final_table.txt
awk '{OFS="\t";print $1,$2,$3,$4,$5}' ${SNPs} > SNPs.txt

Rscript ${script_path}/extract_miRNA.R

bedtools intersect -wo -a pre_miRNA.bed -b SNP.bed | awk 'BEGIN{OFS="\t";print "##fileformat=VCFv4.2""\n""#CHROM","POS","ID","REF","ALT"}{gsub(/T/,"U",$13);gsub(/T/,"U",$14);if($6=="+"){print $4,$8-$2,$10,$13,
$14}else{print $4,$3-$8+1,$10,$13,$14}}' - >SNP.vcf

(head -n 2 SNP.vcf && tail -n +3 SNP.vcf | sort -k 1,1 -k 2,2n ) >SNP_sort.vcf
bcftools view SNP_sort.vcf -Oz -o SNP.vcf.gz
bcftools index -f SNP.vcf.gz

if [ `grep -c "Extended_stem_loop_loc" final_table.txt` -ne '0' ];then
    awk -F"\t" 'NR>1{print ">"$2"\n"$6}' final_table.txt >pre_miRNA.fa
else
    awk -F"\t" 'NR>1{print ">"$1"\n"$3}' final_table.txt >pre_miRNA.fa
fi

python sequenceComplement.py pre_miRNA.bed pre_miRNA.fa
bcftools consensus SNP.vcf.gz --fasta-ref finalPremiRNAsMod.fasta -o finalPremiRNAsMod_SNPs.fasta
python Hairpin_plots.py

pdfunite $(ls -v *.pdf) ${premirnaSNP}
