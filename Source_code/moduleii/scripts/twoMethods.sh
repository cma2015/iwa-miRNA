#!/bin/bash
# two_methods.sh

set -x
script_name="two_methods.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "--curpath   --file path"
  echo "--fasta     --fasta file"
  echo "--txtfile   --alignment file"
  echo "--mergefile --result from miRNATranslate"
  echo "--outcorr   --miRNA precursor"
  echo "--outfile   --results from two methods"
  echo "--version   --Version of script"
  echo "Example: ./$script_name "
  exit 1
}
# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}
main(){
    # Parsing options
    ARGS=`getopt -o vh --long version,help,curpath:,fasta:,txtfile:,mergefile:,expfile:,outcorr:,outfile:,outtardata:,outexp: -n "$0" -- "$@"`
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
            --fasta) fasta=$2; shift 2;;
            --txtfile) txtfile=$2; shift 2;;            
            --mergefile) mergefile=$2; shift 2;;
            --expfile) expfile=$2; shift 2;;
            --outcorr) outcorr=$2; shift 2;;
            --outfile) outfile=$2; shift 2;;
            --outtardata) outtardata=$2; shift 2;;
            --outexp) outexp=$2; shift 2;;
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
    OUTDIR=${curdir}/tmp/222222 #`date +%s%N | cut -c6-13`
    if [ ! -d $OUTDIR ]; then
        mkdir -p ${curdir}/tmp $OUTDIR
    fi

    # python ${SCRIPTPATH}/HTcriteria.py -c ${OUTDIR} -f ${fasta} -t ${txtfile} -m ${mergefile}

    # Rscript ${SCRIPTPATH}/HTcriteria.R ${OUTDIR}/out_pool_merge.txt ${OUTDIR}/pc_criteria.txt

    # cp -r ${SCRIPTPATH}/featureScripts ${OUTDIR}/featureScripts
    cd ${OUTDIR}/featureScripts 
    # python 0extract_features.py
    #Rscript ${SCRIPTPATH}/one-class-svm.R  ${OUTDIR}/00feature_out.txt ${OUTDIR}/out_pool_merge.txt 0.3 ${OUTDIR}/ML_result.txt
    cd ${OUTDIR}
    cp ${SCRIPTPATH}/igraph_rename.sh igraph_rename.sh
    cp ${SCRIPTPATH}/igraph_rename.R igraph_rename.R
    sh igraph_rename.sh
    mkdir rpmData
    tar -zxvf $expfile -C rpmData
    python ${SCRIPTPATH}/expCount.py
    Rscript ${SCRIPTPATH}/expMat.R
    Rscript ${SCRIPTPATH}/mergeOutput.R pc_criteria.txt ML_result.txt out_pool_merge.txt NameChange_out.txt expressionMat.txt final_table.txt
    
    cp final_table.txt ${outfile}
    cp corr_pre-mirna_seq.txt ${outcorr}
    cp expressionMat.txt ${outexp}
    tar -zcvf ${outtardata} -C ${OUTDIR}/data . 
}
main "$@"
