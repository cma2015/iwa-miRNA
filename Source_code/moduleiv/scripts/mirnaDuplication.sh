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
--mirfasta ~/sRNAbox/tools/moduleiv/tmp/mir.gff3 \
--mirbed ~/sRNAbox/tools/moduleiv/tmp/anno.txt \
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
    ARGS=`getopt -o vh --long version,help,curpath:,pcgfasta:,pcggff:,mirfasta:,mirbed:,output1:,output2:,output3: -n "$0" -- "$@"`
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
            --mirfasta) mirfasta=$2; shift 2;;
            --mirtxt) mirtxt=$2; shift 2;;
            --mirgff) mirgff=$2; shift 2;;
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
    # PATH="$PATH:${SCRIPTPATH}"

    # Define the output directory
    OUTDIR=${curdir}/tmp/20200110 #2019_`date +%s%N | cut -c6-13`
    if [ ! -d $OUTDIR ]; then
        mkdir $OUTDIR
    fi

    # genes and proteins
    diamond makedb --in $pcgfasta -d ${OUTDIR}/PCGs
    diamond blastp -d ${OUTDIR}/PCGs -q $pcgfasta -o ${OUTDIR}/matches.m8
    python ${SCRIPTPATH}/GFF3forGP.py $pcggff ${OUTDIR}/GPco.txt ${OUTDIR}/PCGs.gff
    awk 'NR==FNR{a[$1]=$2;next}{$1=a[$1];$2=a[$2];print $0}' ${OUTDIR}/GPco.txt ${OUTDIR}/matches.m8 >${OUTDIR}/PCGs.blast

    # pre-miRNAs
    makeblastdb -in $mirfasta -dbtype nucl -parse_seqids -out ${OUTDIR}/pre-miRNAs
    blastn -query $mirfasta -db ${OUTDIR}/pre-miRNAs -out ${OUTDIR}/ts_mirna.blast -evalue 1 -num_threads 2 -outfmt 6

    ## merge
    cat ${OUTDIR}/PCGs.gff $mirtxt >${OUTDIR}/input.gff
    cat ${OUTDIR}/PCGs.blast ${OUTDIR}/ts_mirna.blast >${OUTDIR}/input.blast

    ${SCRIPTPATH}/MCScanX/MCScanX ${OUTDIR}/input
    ${SCRIPTPATH}/MCScanX/duplicate_gene_classifier ${OUTDIR}/input
}

main("$@")
