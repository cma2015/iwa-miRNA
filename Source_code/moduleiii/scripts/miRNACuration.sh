#!/bin/bash
# two_methods.sh

set -x
script_name="two_methods.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "--curpath     --file path"
  echo "--annoinput   --annotation file"
  echo "--related     --related file"
  echo "--htmlout     --HTML report"
  echo "--version     --Version of script"
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
    ARGS=`getopt -o vh --long version,help,curpath:,annoinput:,related:,speciesInfo:,speciesPath:,htmlout: -n "$0" -- "$@"`
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
            --annoinput) annoinput=$2; shift 2;;
            --related) related=$2; shift 2;;
            --speciesInfo) speciesInfo=$2; shift 2;; 
            --speciesPath) speciesPath=$2; shift 2;; 
            --htmlout) htmlout=$2; shift 2;;
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
    OUTDIR=${curdir}/tmp/000000 #`date +%s%N | cut -c6-13`
    if [ ! -d $OUTDIR ]; then
        mkdir -p ${curdir}/tmp $OUTDIR
    fi

    cd ${OUTDIR} 
    # cp $annoinput final_table.txt
    # cp ${SCRIPTPATH}/miRPage.Rmd miRPage.Rmd
    # cp ${SCRIPTPATH}/miReport.Rmd miReport.Rmd
    # Rscript -e "rmarkdown::render('miRPage.Rmd')"

    # mkdir -p data miRNA_out figure-html
    # tar -zcvf $related -C ${OUTDIR}/data .

    # cut -f 8 final_table.txt | awk 'NR>1' | sed "s/:/-/g" | while read line
    # do
    # echo ${line}

    # sed "s/xxxxxx/${line}/" miReport.Rmd > ${line}.Rmd

    # Rscript -e "rmarkdown::render('${line}.Rmd')"

    # sed -i "s/${line}_files/js_css/" ${line}.html
    # sed -i "s/js_css\/figure-html/figure-html\/${line}/" ${line}.html

    # mkdir figure-html/${line}
    # mv ${line}.html miRNA_out/${line}.html

    # if [ ! -d miRNA_out/js_css ];then
    #     mkdir miRNA_out/js_css
    #     mv -r ${line}_files/* miRNA_out/js_css
    # fi
    # rm ${line}.Rmd
    # rm -r ${line}_files/
    # done

    cp miRPage.html $htmlout
    mkdir -p ${htmlout%.*}_files
    cp -r miRNA_out/ ${htmlout%.*}_files/

}
main "$@"
