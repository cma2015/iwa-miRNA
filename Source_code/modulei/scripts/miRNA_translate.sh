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
  echo "--outhtml    --output html file"
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
    ARGS=`getopt -o vh --long version,help,curpath:,dbfile:,pdfile:,species:,outxt:,outhtml:,outgff: -n "$0" -- "$@"`
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
            --species) species=$2; shift 2;;
            --outxt) outxt=$2; shift 2;;
            --outhtml) outhtml=$2; shift 2;;
            --outgff) outgff=$2; shift 2;;
            -v|--version) version;;
            -h|--help) usage;;
            --) shift; break;;
            *)  echo "Internal error!"; exit 1;;
        esac
    done

    # Add environment path
    script_path=${curdir}/scripts
    PATH="$PATH:${script_path}"

    if [[ -z $dbfile ]] && [[ -z $pdfile ]]; then
        echo "!!! Please check for mandatory options !!!"
        usage
    elif [[ -z $dbfile ]];then
        cp $pdfile $outxt
    else
        ## Define the output directory
        out_dir=${outxt%.*}
        mkdir -p ${out_dir} ${out_dir}/miRNATranslate

        cd ${out_dir}/miRNATranslate
        species_path=${curdir}/`grep -oP "Index/.+" ${species}`
        echo "${species_path}" >speciespath.txt
        ## Data processing
        unzip ${dbfile}
        cp ${curdir}/scripts/gmap_change.py ./gmap_change.py
        if [[ -f ./miRBase.txt ]];then
            bash ${script_path}/version_change.sh --path ${out_dir}/miRNATranslate --genome ${species_path} --flatrna ./miRBase.txt
        fi
        if [[ -f ./PmiREN.txt ]];then
            bash ${script_path}/version_change.sh --path ${out_dir}/miRNATranslate --genome ${species_path} --flatrna ./PmiREN.txt
        fi
        if [[ -f ./sRNAanno.txt ]];then
            bash ${script_path}/version_change.sh --path ${out_dir}/miRNATranslate --genome ${species_path} --flatrna ./sRNAanno.txt
        fi
        if [[ -f ./PlantsmallRNAgenes.txt ]];then
            mv ./PlantsmallRNAgenes.txt ./PlantsmallRNAgenes.txt.bak
            awk '{OFS="\t";print $1,$2,$3,$4,$5,$7,$8,$9,$10,$11,$12}' ./PlantsmallRNAgenes.txt.bak >./PlantsmallRNAgenes.txt
            bash ${script_path}/version_change.sh --path ./ --genome ${species_path} --flatrna ./PlantsmallRNAgenes.txt
        fi
        ## predicted miRNAs
        if [[ ! -z ${pdfile} ]];then
            mkdir ${out_dir}/miRNAPredict
            cp -r ${pdfile%.*}/4rpmData/ ${out_dir}/miRNAPredict/
            cp ${pdfile%.*}/reads* ${out_dir}/miRNAPredict/
            cp ${pdfile%.*}/mapping.bam ${out_dir}/miRNAPredict/
            cp ${pdfile} prediction.txt
        fi

        ## Merge miRNAs
        Rscript ${script_path}/db_merge.R ${out_dir}/miRNATranslate ${species_path}

        ## Detailed report of each miRNAs in html file

        sed "s/<i class='correct'><\/i>/Yes/g;s/<i class='incorrect'><\/i>/No/g" ./Translate_result.txt >./Translate_out.txt
        sed "s/<i class='correct'><\/i>/√/g;s/<i class='incorrect'><\/i>/×/g" ./Translate_result.txt > ./Translate_change.txt

        mkdir -p ./png
        cp ${script_path}/db_page.Rmd ./db_page.Rmd

        python ${script_path}/db_page.py --outpath ${out_dir}/miRNATranslate
        Rscript -e "rmarkdown::render('db_page.Rmd', quiet = T)"

        mkdir -p ${outhtml%.*}_files
        mv ./png ${outhtml%.*}_files/png

        cp ./Translate_change.txt $outxt
        cp db_page.html ${outhtml}
        mv db_page_files ${outhtml%.*}_files/db_page_files
    fi

}
main "$@"
