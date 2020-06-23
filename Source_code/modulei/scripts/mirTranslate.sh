#!/bin/bash
# mirTranslate.sh

set -x
script_name="mirTranslate.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "--curpath   --file path"
  echo "--dbfile    --miRNA candidates from databases"
  echo "--pdfile    --predicted miRNAs"
  echo "--outxt     --output file"
  echo "--outgff    --output GFF3 file"  
  echo "--version   --Version of script"
  echo "Example: sh ~/sRNAbox/tools/modulei/scripts/mirTranslate.sh --curpath ~/sRNAbox/tools/modulei \
  --dbfile ~/sRNAbox/tools/modulei/tmp/QgeCW9Yyuo/00merge.txt \
  --pdfile ~/sRNAbox/tools/modulei/tmp/20200206/out_info.txt \
  --outxt ~/sRNAbox/tools/modulei/tmp/db_pd_merge.txt \
  --outgff ~/sRNAbox/tools/modulei/tmp/db_pd_merge.gff3"
  exit 1
}
# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}
main(){
    # Parsing options
    ARGS=`getopt -o vh --long version,help,curpath:,dbfile:,pdfile:,outxt:,outgff: -n "$0" -- "$@"`
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
            --dbfile) dbfile=$2; shift 2;;
            --pdfile) pdfile=$2; shift 2;;
            --outxt) outxt=$2; shift 2;;
            --outgff) outgff=$2; shift 2;;
            -v|--version) version;;
            -h|--help) usage;;
            --) shift; break;;
            *)  echo "Internal error!"; exit 1;;
        esac
    done

    # Add environment path
    SCRIPTPATH=${curdir}/scripts
    PATH="$PATH:${SCRIPTPATH}"

    # Define the output directory
    OUTDIR=${curdir}/tmp/111222 #`date +%s%N | cut -c6-13`
    if [ ! -d $OUTDIR ]; then
        mkdir $OUTDIR
    fi

    Rscript ${SCRIPTPATH}/mir_translate.R $pdfile $dbfile ${OUTDIR}/db_pd_merge.txt
    cp ${OUTDIR}/db_pd_merge.txt $outxt

}
main "$@"
