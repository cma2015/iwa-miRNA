#!/bin/bash
# mirnaAssociation.sh

set -x
script_name="mirnaAssociation.sh"
script_ver="1.0.0"

# Help function
usage() {
  echo "--curpath     --file path"
  echo "--miranno     --miRNA annotation"
  echo "--genomeanno  --Annotation of genome"
  echo "--rmd         --Rmarkdown report"
  echo "--version     --Version of script"
  echo "Example: bash ~/sRNAbox/tools/moduleiv/scripts/mirnaAssociation.sh --path ~/sRNAbox/tools/moduleiv \
--mirtxt ~/sRNAbox/tools/moduleiv/tmp/mir.gff3 --genomeanno ~/sRNAbox/tools/moduleiv/tmp/anno.txt \
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
    OUTDIR=${curdir}/../tmp/associationAnalysis_`date +%s%N | cut -c6-13`
    if [ ! -d $OUTDIR ]; then
        mkdir -p ${curdir}/tmp $OUTDIR
    fi
    cd ${OUTDIR}
    Rscript ${SCRIPTPATH}/mirExtract.R ${mirtxt} ${OUTDIR}/mir.txt
    awk '{OFS="\t";print $1,$2,$3,$4,".",$5,$6,$7}' ${genomeanno} >genome_anno.txt
    bedtools intersect -wo -a ${OUTDIR}/mir.txt -b genome_anno.txt >${OUTDIR}/miRNA_intersect_out.txt
    cp ${SCRIPTPATH}/mirAssociation.Rmd ${OUTDIR}/mirAssociation.Rmd
    Rscript -e "rmarkdown::render('mirAssociation.Rmd')"
    cp ${OUTDIR}/mirAssociation.html ${rmd}
}

main $@
