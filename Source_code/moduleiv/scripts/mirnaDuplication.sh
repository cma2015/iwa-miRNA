#!/bin/bash
# mirnaDuplication.sh

script_name="mirnaDuplication.sh"
script_ver="1.0.0"

# Help function
usage() {
  echo "--curpath     --file path"
  echo "--pcgfasta    --miRNA GFF3"
  echo "--pcggff      --Annotation of genome"
  echo "--mirfasta    --mirfasta"
  echo "--mirbed      --mirbed"
  echo "--output1     --output1"
  echo "--output2     --output2"
  echo "--output3     --output3"
  echo "--version     --version of script"
  echo "Example: bash ~/sRNAbox/tools/moduleiv/scripts/mirnaDuplication.sh \
--curpath ~/sRNAbox/tools/moduleiv \
--pcgfasta ~/sRNAbox/tools/moduleiv/tmp/mir.gff3 \
--pcggff ~/sRNAbox/tools/moduleiv/tmp/anno.txt \
--mirAnno ~/sRNAbox/tools/moduleiv/tmp/final_table.txt \
--output1 ~/sRNAbox/tools/moduleiv/tmp/mir.gff3 \
--output2 ~/sRNAbox/tools/moduleiv/tmp/anno.txt \
--output3 out_test.html"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){
    ## Parsing options
    ARGS=`getopt -o vh --long version,help,curpath:,pcgfasta:,pcggff:,mirAnno:,output1:,output2:,output3: -n "$0" -- "$@"`
    if [ $? != 0 ]; then
        echo "Terminating..."
        exit 1
    fi
    eval set -- "${ARGS}"
    echo formatted parameters=[$@]
    while true
    do
        case "$1" in
            --curpath) curdir=$2; shift 2;;
            --pcgfasta) pcgfasta=$2; shift 2;;
            --pcggff) pcggff=$2; shift 2;;
            --mirAnno) mirAnno=$2; shift 2;;
            --output1) output1=$2; shift 2;;
            --output2) output2=$2; shift 2;;
            --output3) output3=$2; shift 2;;
            -v|--version) version;;
            -h|--help) usage;;
            --) shift; break;;
            *)  echo "Internal error!"; exit 1;;
        esac
    done

    ## Add environment path
    SCRIPTPATH=${curdir}/scripts

    # Define the output directory
    OUTDIR=${curdir}/../tmp/duplicationMIR_`date +%s%N | cut -c6-13`

    if [ ! -d $OUTDIR ]; then
        mkdir $OUTDIR
    fi

    cd ${OUTDIR}
    # genes and proteins
    diamond makedb --in $pcgfasta -d PCGs
    diamond blastp -d PCGs -q $pcgfasta -o matches.m8
    python ${SCRIPTPATH}/GFF3forGP.py $pcggff GPco.txt PCGs.gff
    awk 'NR==FNR{a[$1]=$2;next}{OFS="\t";$1=a[$1];$2=a[$2];print $0}' GPco.txt matches.m8 >PCGs.blast

    # pre-miRNAs
    cp ${mirAnno} final_table.txt

    Rscript ${SCRIPTPATH}/mirFastaTxt.R final_table.txt mir.gff mir.fasta

    makeblastdb -in mir.fasta -dbtype nucl -parse_seqids -out pre-miRNAs
    blastn -query mir.fasta -db pre-miRNAs -out ts_mirna.blast -evalue 1 -num_threads 2 -outfmt 6

    ## merge
    cat PCGs.gff mir.gff >input.gff
    cat PCGs.blast ts_mirna.blast >input.blast

    ${SCRIPTPATH}/MCScanX/MCScanX input
    ${SCRIPTPATH}/MCScanX/duplicate_gene_classifier input

    Rscript ${SCRIPTPATH}/R_syntenic_regions.R

    cp summary_of_miRNAs_duplication.txt ${output1}
    cp Genomic_duplication_and_miRNA_expansion.pdf ${output2}
    cp miRNA_duplication_in_blocks.txt ${output3}
}

main "$@"
