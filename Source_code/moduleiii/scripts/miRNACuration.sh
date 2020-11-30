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
    ARGS=`getopt -o vh --long version,help,curpath:,annoinput:,related:,sampleInfo:,description:,genomeanno:,matexp:,htmlout: -n "$0" -- "$@"`
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
            --sampleInfo) sampleInfo=$2; shift 2;;
            --description) description=$2; shift 2;;
            --genomeanno) genomeanno=$2; shift 2;;
            --htmlout) htmlout=$2; shift 2;;
            --matexp) matexp=$2; shift 2;;
            -v|--version) version;;
            -h|--help) usage;;
            --) shift; break;;
            *)  echo "Internal error!"; exit 1;;
        esac
    done

    # Add environment path
    script_path=${curdir}/scripts
    PATH="$PATH:${script_path}"

    # Define the output directory
    out_dir=${annoinput%.*}
    if [ ! -d $out_dir/manualCuration ]; then
        mkdir -p $out_dir/manualCuration
    fi

    cd ${out_dir}/manualCuration
    mkdir miRNA_out

    cp ${annoinput} final_table.txt
    awk '{OFS="\t";print $1,$2}' ${sampleInfo} > sample_info.txt
    cp ${out_dir}/miRNASelection/miRNA_in_sample.txt miRNA_in_sample.txt
    if [ -f $description ];then
        cp ${description} gene_description.txt
    fi
    cp ${script_path}/miRPage.Rmd miRPage.Rmd
    cp ${script_path}/miReport.Rmd miReport.Rmd
    
    ## genomic sources
    if [[ -f ${genomeanno} ]];then
        Rscript ${script_path}/miRNA_extraction.R final_table.txt miRNA_list.txt
        awk '{OFS="\t";print $1,$2,$3,$4,".",$5,$6,$7}' ${genomeanno} >genome_anno.txt
        awk '$7~/[tT][eE]/' genome_anno.txt > TEs.txt
        awk '$7!~/[tT][eE]/' genome_anno.txt > non_TEs.txt
        bedtools intersect -wo -a miRNA_list.txt -b TEs.txt -f 0.5 > miRNA_TE.txt
        bedtools intersect -wo -a miRNA_list.txt -b non_TEs.txt -f 0.5 -s > miRNA_non_TE.txt
        Rscript ${script_path}/genomic_source.R
    else
        Rscript ${script_path}/genomic_no_source.R
    fi

    ##
    Rscript -e "rmarkdown::render('miRPage.Rmd', quiet=T)"

    awk -F"\t" 'NR>1{print $5}' final_table.txt | sed "s/:/-/g" >tmp_name.txt
    awk -F"\t" 'NR>1{print ">"$5"_5p\n"$11"\n>"$5"_3p\n"$14}' final_table.txt >mature_miRNAs.fa

    samtools sort ../miRNAPredict/mapping.bam -o mapping.bam
    samtools index mapping.bam

    speciesgenome=`cat ../miRNATranslate/speciespath.txt`/Genome/Genome.fa
    python ${script_path}/struVis_extraction.py ${script_path} mapping.bam ${speciesgenome} final_table.txt

    # speciespath=`cat ${out_dir}/miRNATranslate/speciespath.txt`

    #if [ -s psRNAtarget_MIT.out ];then
    #    if [ `grep -c "miRNA_Acc." psRNAtarget_MIT.out` -ne '0' ];then
    #        echo "pass!"
    #    else
    #        python ${script_path}/PsRNAtarget.py ${speciespath} ps_path.sh
    #    fi
    #else
    #    python ${script_path}/PsRNAtarget.py ${speciespath} ps_path.sh
    #fi

    if [ -s ../miRNASelection/ps_path.sh ];then
        sh ../miRNASelection/ps_path.sh
    fi

    # awk -F"\t" 'NR>1{print $5}' final_table.txt | while read line
    # do
    #     grep "${line}" psRNAtarget_MIT.out >${out_dir}/miRNASelection/data/${line}.mti
    # done
    awk -v RS=">" 'NR>1{print $1"\t"$2}' ../miRNASelection/mature_miRNAs.fa >mature_miRNAs.txt
    Rscript ${script_path}/PsRNAtarget_extraction.R

    cat tmp_name.txt | while read line
    do
        echo ${line}
        sed "s/xxxxxx/${line}/" miReport.Rmd > ${line}.Rmd
        Rscript -e "rmarkdown::render('${line}.Rmd', quiet=TRUE)"
        sed -i "s/${line}_files/js_css/" ${line}.html

        mv ${line}.html miRNA_out/${line}.html
        if [ ! -d miRNA_out/js_css ];then
            mkdir miRNA_out/js_css
	    cp -r ${line}_files/* miRNA_out/js_css
        else
            rm -r ${line}_files
        fi
        rm ${line}.Rmd
    done

    mv miRPage.html $htmlout
    mkdir -p ${htmlout%.*}_files
    mv miRPage_files ${htmlout%.*}_files
    mv miRNA_out ${htmlout%.*}_files
    
    runtime=$(date "+%Y-%m-%d %H:%M:%S")
echo """

${runtime}: ${out_dir}

""" >> /home/runlog.txt

}

main "$@"
