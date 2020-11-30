#!/bin/bash
# miRNA_association.sh

set -x
script_name="miRNA_association.sh"
script_ver="1.0.0"

# Help function
usage() {
  echo "--curpath     --file path"
  echo "--miranno     --miRNA annotation"
  echo "--genomeanno  --Annotation of genome"
  echo "--rmd         --Rmarkdown report"
  echo "--version     --Version of script"
  echo "Example: bash ./moduleiv/scripts/miRNA_association.sh --path ./tools/moduleiv \
--mirtxt ./tools/moduleiv/tmp/mir.gff3 --genomeanno ./tools/moduleiv/tmp/anno.txt \
--rmd out_test.html"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){
    ## Parsing options
    ARGS=`getopt -o vh --long version,help,curpath:,miranno:,genomeanno:,rmd: -n "$0" -- "$@"`
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
            --miranno) mirtxt=$2; shift 2;;
            --genomeanno) genomeanno=$2; shift 2;;
            --filter) filter=$2; shift 2;;
            --rmd) rmd=$2; shift 2;;
            -v|--version) version;;
            -h|--help) usage;;
            --) shift; break;;
            *)  echo "Internal error!"; exit 1;;
        esac
    done

    ## Add environment path
    SCRIPTPATH=${curdir}/scripts
    PATH="$PATH:${SCRIPTPATH}"

    # Define the output directory
    OUTDIR=${mirtxt%.*}/associationAnalysis

    if [ ! -d ${mirtxt%.*} ]; then
        mkdir ${mirtxt%.*}
    fi

    if [ ! -d $OUTDIR ]; then
        mkdir -p $OUTDIR
    fi
    
    cd ${OUTDIR}
    Rscript ${SCRIPTPATH}/miRNA_extraction.R ${mirtxt} miRNA_list.txt
    awk '{OFS="\t";print $1,$2,$3,$4,".",$5,$6,$7}' ${genomeanno} >genome_anno.txt
    awk '$7~/[tT][eE]/' genome_anno.txt > TEs.txt
    awk '$7!~/[tT][eE]/' genome_anno.txt > non_TEs.txt
    bedtools intersect -wo -a miRNA_list.txt -b TEs.txt -f 0.5 > miRNA_TE.txt
    bedtools intersect -wo -a miRNA_list.txt -b non_TEs.txt -f 0.5 -s > miRNA_non_TE.txt
    # bedtools intersect -wao -a miRNA_list.txt -b genome_anno.txt -f 0.5 >miRNA_intersect_out.txt
    cp ${SCRIPTPATH}/miRNA_genome_association.Rmd miRNA_genome_association.Rmd
    Rscript -e "rmarkdown::render('miRNA_genome_association.Rmd', quiet = T)"
    mkdir -p ${rmd%.*}_files
    mv miRNA_genome_association_files ${rmd%.*}_files/miRNA_genome_association_files
    cp miRNA_genome_association.html ${rmd}
}

main $@
