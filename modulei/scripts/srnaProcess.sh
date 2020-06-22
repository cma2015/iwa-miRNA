#!/bin/bash

set -x
script_name="srnaProcess.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  echo "-f  --The current path."
  echo "-s  --SRA accession numbers."
  echo "-u  --User's data"
  echo "-a  --Adatpter file"
  echo "-q  --Quality score"
  echo "-m  --Minimal read length"
  echo "-n  --Maximum read length"
  echo "-c  --Minimal read number for clean files"
  echo "-g  --tar.gz file"
  echo "-i  --multiQC report"
  echo "-l  --Rmarkdown report"
  echo "-p  --Number of threads"
  echo "-v  --Version of script"
  echo "Example: ./$script_name -s 'SRR1039523,SRR1039524' -u 'maize_samples.tar.gz' -a 'adpter.txt' -q '20' -m '18' -n '26' -c '1000000' -p '3'"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){
    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:s:u:a:q:m:n:c:p:g:i:l:hv opt
        do
            case $opt in
                f) curdir=$OPTARG;;
                s) sralist=$OPTARG;;
                u) userdata=$OPTARG;;
                a) adatper=$OPTARG;;
                q) quality=$OPTARG;;
                m) min=$OPTARG;;
                n) max=$OPTARG;;
                c) readnum=$OPTARG;;
                p) thread=$OPTARG;;
                g) outputgz=$OPTARG;;
                i) multiQCres=$OPTARG;;
                l) outputhtml=$OPTARG;;
                h) usage;;
                v) version;;
            esac
        done
    shift $(($OPTIND -1))

    ## Add environment path
    scriptpath=${curdir}/scripts
    PATH="$PATH:${scriptpath}"

    # Check for mandatory options
    if [[ -z $quality ]] || [[ -z $min ]] || [[ -z $max ]]; then
        echo "!!! Please check for mandatory options !!!"
        usage
    fi

    # Define the output directory
    out_dir=${curdir}/tmp/000111 #2019_`date +%s%N | cut -c6-13`

    if [ ! -d $out_dir ]; then
        mkdir -p $out_dir $out_dir/00rawdata
    fi

    # Write sample name to yaml
    if [ ! -z $userdata ]; then
        tar xvf ${userData} -C ${out_dir}/00rawdata
        for folder in ${out_dir}/00rawdata/*
        do
            if [ -d $folder ];then
                for file in ${folder}/*
                do
                    filename=`basename ${file}`
                    mv ${file} ${out_dir}/00rawdata/${filename}
                done
                rmdir $folder
            fi
        done
    fi

    if [ ! -z $sralist ]; then
        array_sra=(${sralist//[__cn__|,]/ })
        echo -e "inputpath:\n    ${curdir}:" >${out_dir}/config.yml
        echo -e "outpath:\n    ${out_dir}/00rawdata:\nsamples:" >>${out_dir}/config.yml
        for easra in ${array_sra[@]}
            do
                echo "    ${easra}:" >>${out_dir}/config.yml
            done
        cat ${out_dir}/config.yml
        snakemake --snakefile ${curdir}/scripts/srnaDownload.smk --configfile ${out_dir}/config.yml --cores ${thread} --latency-wait 120
    fi

    if [ -f ${out_dir}/Adapter.txt ];then
        rm ${out_dir}/Adapter.txt
    fi
    if [ ! -z $adatper ]; then
        cat ${adatper} | while read line
        do
        srr=`echo $line | cut -d " " -f 1`
        ada=`echo $line | cut -d " " -f 2`
        phred=`zcat ${out_dir}/00rawdata/${srr}.fastq.gz | fastq_phred.pl - | head -n 1`
        barcode=`echo $line | awk '{if($3~/[ATCG]/){print length($3)+1}else{print $3}}'`
        echo -e "${srr}\t$ada\t$phred\t$barcode" >>${out_dir}/Adapter.txt
        done
    else
        for file in ${out_dir}/00rawdata/*.fastq.gz
        do
            phred=`zcat ${file} | fastq_phred.pl - | head -n 1`
            barcode=`zcat ${file} | awk '{if((NR+2)%4==0){a[substr($1,1,2)]++;b[substr($1,1,3)]++;c[substr($1,1,4)]++}}\
            END{if(length(c)==1){print 5}else if(length(b)==1){print 4}else \
			if(length(a)==1){print 3}else{print "-"}}' - `
            file_name=`basename ${file}`
            if [ -f ${out_dir}/tmp_seq.txt ];then
                rm ${out_dir}/tmp_seq.txt
            fi
            zcat ${file} | grep "TGACAGAAGAGAGTGAGCAC" - >${out_dir}/tmp_seq.txt
            adapter_py=`dataAdapter.py ${out_dir}/tmp_seq.txt`
            if [ -z $adapter_py ];then adapter_py='NA' ;fi
            adapter_dnapi=`zcat ${file} | dnapi.py -`
            if [ $adapter_py != 'NA' ];then
                echo -e "${file_name%%.*}\t$adapter_py\t$phred\t$barcode" >>${out_dir}/Adapter.txt
            else
                echo -e "${file_name%%.*}\t$adapter_dnapi\t$phred\t$barcode" >>${out_dir}/Adapter.txt
            fi
        done
    fi
    
    # Data processing
    echo -e "inputpath:\n    $curdir:" >${out_dir}/config.yml
    echo -e "quality:\n    '$quality':" >>${out_dir}/config.yml
    echo -e "min:\n    '$min':" >>${out_dir}/config.yml
    echo -e "max:\n    '$max':" >>${out_dir}/config.yml
    echo -e "outpath:\n    $out_dir:\nsamples:" >>${out_dir}/config.yml
    ls $out_dir/00rawdata/*.fastq.gz | sed -r 's/.+\///;s/\..+//' | while read line
    do
        echo "    ${line}:" >>${out_dir}/config.yml
    done
    snakemake --snakefile ${curdir}/scripts/srnaProcess.smk --configfile ${out_dir}/config.yml --cores ${thread} --latency-wait 120
    
    multiqc -q -o ${out_dir}/1multiQC -d ${out_dir}/1fastqc -dd 1 -f -v -n multiQC -i sRNA_seqQC
    # Show the results
    awk -v readnum=${readnum} '$5<readnum' ${out_dir}/00Table_Summary_of_sRNA-seq_data.txt | cut -f 1 | while read line
    do
        rm ${out_dir}/2collapsedata/${line}.fasta ${out_dir}/3readlength/${line}_all.txt ${out_dir}/3readlength/${line}_unique.txt
    done
    cat ${out_dir}/3readlength/*.txt > ${out_dir}/Length_count.txt
    
    tar -zcvf $outputgz -C ${out_dir}/2collapsedata .
    cp ${out_dir}/1multiQC/multiQC.html ${multiQCres}
    
    # Rmarkdown report
    cp ${curdir}/scripts/srnaProcess.Rmd ${out_dir}
    Rscript -e "rmarkdown::render('${out_dir}/srnaProcess.Rmd')"
    cp ${out_dir}/srnaProcess.html ${outputhtml}
}

main "$@"
