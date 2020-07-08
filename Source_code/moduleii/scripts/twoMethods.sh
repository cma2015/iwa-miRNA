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
  echo "--mirmerge --result from miRNATranslate"
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
    ARGS=`getopt -o vh --long version,help,curpath:,mirmerge:,expfile:,stemloop:,structure:,abias:,sbias:,minlen:,maxlen:,positive:,nuval:,outfile: -n "$0" -- "$@"`
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
            --mirmerge) mirmerge=$2; shift 2;;
            --expfile) expfile=$2; shift 2;;
            --stemloop) stemloop=$2; shift 2;;
            --structure) structure=$2; shift 2;;
            --abias) abias=$2; shift 2;;
            --sbias) sbias=$2; shift 2;;
            --minlen) minlen=$2; shift 2;;
            --maxlen) maxlen=$2; shift 2;;
            --positive) positive=$2; shift 2;;
            --nuval) nuval=$2; shift 2;;
            --outfile) outfile=$2; shift 2;;
            -v|--version) version;;
            -h|--help) usage;;
            --) shift; break;;
            *)  echo "Internal error!"; exit 1;;
        esac
    done

    # Add environment path
    script_path=${curdir}/scripts

    # Define the output directory
    out_dir=${outfile%.*}

    mkdir -p ${out_dir} ${out_dir}/miRNASelection

    cp -r ${mirmerge%.*}/* ${out_dir}/

    cd ${out_dir}/miRNASelection

    python ${script_path}/HTcriteria.py -c ${out_dir}/miRNASelection -f ../miRNAPredict/reads_18_26.fa -t ../miRNAPredict/reads_18_26.txt -m ${mirmerge}
    Rscript ${script_path}/HTcriteria.R ${out_dir}/miRNASelection/out_pool_merge.txt ${out_dir}/miRNASelection/pc_criteria.txt ${stemloop} ${structure} ${abias} ${sbias} ${minlen} ${maxlen}

    cp -r ${script_path}/featureScripts ${out_dir}/miRNASelection/featureScripts
    cd ${out_dir}/miRNASelection/featureScripts
    python 0extract_features.py

    cd ${out_dir}/miRNASelection
    if [ ! -z ${positive} ];then
        cp ${positive} positive_name.txt
        Rscript ${script_path}/one-class-svm.R  00feature_out.txt out_pool_merge.txt ${nuval} ML_result.txt positive_name.txt
    else
        Rscript ${script_path}/one-class-svm.R  00feature_out.txt out_pool_merge.txt ${nuval} ML_result.txt
    fi
    cp ${script_path}/igraph_rename.sh igraph_rename.sh
    cp ${script_path}/igraph_rename.R igraph_rename.R
    sh igraph_rename.sh

    python ${script_path}/expCount.py ${out_dir}/miRNAPredict
    Rscript ${script_path}/expMat.R
    Rscript ${script_path}/mergeOutput.R pc_criteria.txt ML_result.txt out_pool_merge.txt NameChange_out.txt expressionMat.txt final_table.txt

    awk -F"\t" 'NR>1{print $11"\n"$14}' final_table.txt | sort -u >mature_miRNAs.fa
    python ${script_path}/matureMat.py ${out_dir}/miRNAPredict
    Rscript ${script_path}/matureMat.R

    ## Output
    cd ${out_dir}/miRNASelection
    cp final_table.txt ${outfile}

    # cp corr_pre-mirna_seq.txt ${outcorr}
    # cp expressionMat.txt ${outexp}
    # tar -zcvf ${outtardata} -C ./data/ .
    # cp miRNA_in_sample.txt ${matexp}
}
main "$@"
