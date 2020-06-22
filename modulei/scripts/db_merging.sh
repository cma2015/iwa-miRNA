#!/bin/bash
# db_merging.sh

set -x
script_name="db_merging.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "--curpath   --file path"
  echo "--rnatype   --rna type"
  echo "--species   --species"
  echo "--mirbase   --use mirbase"
  echo "--mirbaseC  --change id in mirbase"
  echo "--PmiREN    --use PmiREN"
  echo "--PmiRENC   --change id in PmiREN"
  echo "--sRNAanno  --use sRNAanno"
  echo "--sRNAannoC --change id in RPM"
  echo "--Pgenes    --use PlantsmallRNAgenesS"
  echo "--PgenesC   --change id in PlantsmallRNAgenesS"
  echo "--genome    --genome path"
  echo "--outxt     --merged file"
  echo "--rmd       --Rmarkdown report"
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
    ARGS=`getopt -o vh --long version,help,curpath:,rnatype:,species:,mirbase:,\
mirbaseC:,PmiREN:,PmiRENC:,sRNAanno:,sRNAannoC:,Pgenes:,PgenesC:,genome:,outxt:,rmd: -n "$0" -- "$@"`
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
            --rnatype) rnatype=$2; shift 2;;
            --species) species=$2; shift 2;;
            --mirbase) mirbase=$2; shift 2;;
            --mirbaseC) mirbaseC=$2; shift 2;;
            --PmiREN) PmiREN=$2; shift 2;;
            --PmiRENC) PmiRENC=$2; shift 2;;
            --sRNAanno) sRNAanno=$2; shift 2;;
            --sRNAannoC) sRNAannoC=$2; shift 2;;
            --Pgenes) Pgenes=$2; shift 2;;
            --PgenesC) PgenesC=$2; shift 2;;
            --genome) genome=$2; shift 2;;
            --outxt) outxt=$2; shift 2;;
            --rmd) rmd=$2; shift 2;;
            -v|--version) version;;
            -h|--help) usage;;
            --) shift; break;;
            *)  echo "Internal error!"; exit 1;;
        esac
    done

    # Add environment path
    SCRIPTDIR=${curdir}/scripts

    dbname=""
    if [[ ! -z $mirbase ]];then
        dbname=miRBase
    fi
    if [[ ! -z $PmiREN ]];then
        if [[ $dbname == "" ]];then
            dbname=PmiREN
        else
            dbname=${dbname},PmiREN
        fi
    fi
    if [[ ! -z $sRNAanno ]];then
        if [[ $dbname == "" ]];then
            dbname=sRNAanno
        else
            dbname=${dbname},sRNAanno
        fi
    fi
    if [[ ! -z $Pgenes ]];then
        if [[ $dbname == "" ]];then
            dbname=PlantsmallRNAgenes
        else
            dbname=${dbname},PlantsmallRNAgenes
        fi
    fi

    ## Extract miRNAs in databases
    out_dir=`python ${SCRIPTDIR}/db_mining.py \
    --curpath ${curdir} --species ${species} --database ${dbname} --outformat txt | grep "modulei"`

    genomepath=`cat ${genome}`

    cp ${curdir}/scripts/gmap_change.py ${out_dir}
    if [[ ! -z $mirbaseC ]];then
        bash ${SCRIPTDIR}/version_change.sh --path ${out_dir} --genome ${curdir}/${genomepath} --flatrna ${out_dir}/miRBase.txt
    fi
    if [[ ! -z $PmiRENC ]];then
        bash ${SCRIPTDIR}/version_change.sh --path ${out_dir} --genome ${curdir}/${genomepath} --flatrna ${out_dir}/PmiREN.txt
    fi
    if [[ ! -z $sRNAannoC ]];then
        bash ${SCRIPTDIR}/version_change.sh --path ${out_dir} --genome ${curdir}/${genomepath} --flatrna ${out_dir}/sRNAanno.txt
    fi
    if [[ ! -z $PgenesC ]];then
        bash ${SCRIPTDIR}/version_change.sh --path ${out_dir} --genome ${curdir}/${genomepath} --flatrna ${out_dir}/PlantsmallRNAgenes.txt
    fi

    ## Merge miRNAs
    Rscript ${SCRIPTDIR}/db_merge.R ${out_dir} ${dbname} ${curdir}/${genomepath}
    sed "s/<i class='correct'><\/i>/√/g;s/<i class='incorrect'><\/i>/×/g" ${out_dir}/00merge.txt > $outxt

    ## Detailed report of each miRNAs in html file
    mkdir -p ${out_dir}/png
    cp ${SCRIPTDIR}/db_page.Rmd ${out_dir}/db_page.Rmd
    Rscript -e "rmarkdown::render('${out_dir}/db_page.Rmd',output_dir='${out_dir}')"
    python ${SCRIPTDIR}/db_page.py --outpath ${out_dir}

    ## Final report file
    cp ${out_dir}/db_page.html ${rmd}
    mkdir -p ${rmd%.*}_files
    cp -r ${out_dir}/png ${rmd%.*}_files/

}
main "$@"
