#!/bin/bash

set -x
script_name="srnaProcess.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){
    # Parsing options
    ARGS=`getopt -o vh --long version,help,curpath:,species:,sralist:,userdata:,adatper:,quality:,min:,max:,\
method:,mirreadnum:,mirmin:,mirmax:,mirmap:,mirsum:,mirleast:,mirsamnum:,thread:,multiqcres:,outputhtml:,outxt:,exptar: -n "$0" -- "$@"`
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
            --species) species=$2; shift 2;;
            --sralist) sralist=$2; shift 2;;
            --userdata) userdata=$2; shift 2;;
            --adatper) adatper=$2; shift 2;;
            --quality) quality=$2; shift 2;;
            --min) min=$2; shift 2;;
            --max) max=$2; shift 2;;
            --method) method=$2; shift 2;;
            --mirreadnum) mirreadnum=$2; shift 2;;
            --mirmin) mirmin=$2; shift 2;;
            --mirmax) mirmax=$2; shift 2;;
            --mirmap) mirmap=$2; shift 2;;
            --mirsum) mirsum=$2; shift 2;;
            --mirleast) mirleast=$2; shift 2;;
            --mirsamnum) mirsamnum=$2; shift 2;;
            --thread) thread=$2; shift 2;;
            --multiqcres) multiqcres=$2; shift 2;;
            --outputhtml) outputhtml=$2; shift 2;;
            --outxt) outxt=$2; shift 2;;
            --exptar) exptar=$2; shift 2;;
            -v|--version) version;;
            -h|--help) usage;;
            --) shift; break;;
            *)  echo "Internal error!"; exit 1;;
        esac
    done


    ## Add environment path
    script_path=${curdir}/scripts
    PATH="$PATH:${script_path}"

    # Check for mandatory options
    if [[ -z $sralist ]] && [[ -z $userdata ]]; then
        echo "!!! Please check for mandatory options !!!"
        usage
    fi

    # Define the output directory
    out_dir=${outxt%.*}
    if [ ! -d $out_dir ]; then
        mkdir -p $out_dir
    fi
    cd ${out_dir} && mkdir -p 00rawdata 4mircatOut 5mirdpOut
    tmp_species=`grep -oP "Index/.+" $species`
    species_path=${curdir}/${tmp_species}

    ## Write sample name to yaml
    if [ ! -z $userdata ] && [ $userdata != "None" ]; then
        cd 00rawdata && unzip ${userdata} && cd ..
        gunzip -k ./00rawdata/*
        # for folder in 00rawdata/*
        # do
        #     if [ -d $folder ];then
        #         for file in ${folder}/*
        #         do
        #             filename=`basename ${file}`
        #             mv ${file} 00rawdata/${filename}
        #         done
        #         rmdir $folder
        #     fi
        # done
    fi

    if [ $thread -gt 5 ];then
        limited=5;
    else
        limited=$thread
    fi

    if [ ! -z "$sralist" ]; then
        sralist=${sralist//./}
        sralist=${sralist// /}
        array_sra=(${sralist//[__cn__|,]/ })
        echo -e "inputpath:\n    ${curdir}:" >config.yml
        echo -e "outpath:\n    00rawdata:\nsamples:" >>config.yml
        for easra in ${array_sra[@]}
            do
                echo "    ${easra}:" >>config.yml
            done
        cat config.yml
        snakemake --snakefile ${curdir}/scripts/srnaDownload.smk --configfile config.yml --cores ${limited} --latency-wait 120
    fi

    if [ -f Adapter.txt ];then
        rm Adapter.txt
    fi

    if [ ! -z $adatper ]; then
        awk -v RS=">" 'NR>1{print $1,$2}' ${adatper} | while read line
        do
        srr=`echo $line | cut -d " " -f 1`
        ada=`echo $line | cut -d " " -f 2`
        phred=`zcat 00rawdata/${srr}.fastq.gz | fastq_phred.pl - | head -n 1`
        barcode=`echo $line | awk '{if($3~/[ATCG]/){print length($3)+1}else{print $3}}'`
        echo -e "${srr}\t$ada\t$phred\t$barcode" >>Adapter.txt
        done
    else
        for file in 00rawdata/*.gz
        do
            phred=`zcat ${file} | fastq_phred.pl - | head -n 1`
            barcode=`zcat ${file} | awk '{if((NR+2)%4==0){a[substr($1,1,2)]++;b[substr($1,1,3)]++;c[substr($1,1,4)]++}}\
            END{if(length(c)==1){print 5}else if(length(b)==1){print 4}else \
			if(length(a)==1){print 3}else{print "-"}}' - `
            file_name=`basename ${file}`
            if [ -f tmp_seq.txt ];then
                rm tmp_seq.txt
            fi
            zcat ${file} | grep "TGACAGAAGAGAGTGAGCAC" - >tmp_seq.txt
            adapter_py=`dataAdapter.py tmp_seq.txt`
            if [ -z $adapter_py ];then adapter_py='NA' ;fi
            adapter_dnapi=`zcat ${file} | dnapi.py -`
            if [ $adapter_py != 'NA' ];then
                echo -e "${file_name%%.*}\t$adapter_py\t$phred\t$barcode" >>Adapter.txt
            else
                echo -e "${file_name%%.*}\t$adapter_dnapi\t$phred\t$barcode" >>Adapter.txt
            fi
        done
    fi

    ## Data processing
    echo -e "inputpath:\n    $curdir:" >config.yml
    echo -e "quality:\n    '$quality':" >>config.yml
    echo -e "min:\n    '$min':" >>config.yml
    echo -e "max:\n    '$max':" >>config.yml
    echo -e "outpath:\n    $out_dir:\nsamples:" >>config.yml
    ls $out_dir/00rawdata/*.gz | sed -r 's/.+\///;s/\..+//' | while read line
    do
        echo "    ${line}:" >>config.yml
    done
    snakemake --snakefile ${curdir}/scripts/srnaProcess.smk --configfile config.yml --cores ${limited} --latency-wait 120
    multiqc -q -o 1multiQC -d 1fastqc -dd 1 -f -v -n multiQC -i sRNA_seqQC

    ## Data processing
    echo -e "inputpath:\n    $curdir:" >./config.yml
    echo -e "speciespath:\n    $species_path:" >>./config.yml
    echo -e "readnum:\n    '$mirreadnum:'" >>./config.yml
    echo -e "multimap:\n    '$mirmap':" >>./config.yml
    echo -e "sum:\n    '$mirsum':" >>./config.yml
    echo -e "least:\n    '$mirleast':" >>./config.yml
    echo -e "samnum:\n    '$mirsamnum':" >>./config.yml
    echo -e "outpath:\n    ${out_dir}:\nsamples:" >>./config.yml

    ls 2collapsedata/ | while read line
    do
        echo "    ${line%.*}:" >>./config.yml
    done

    snakemake --snakefile ${script_path}/mirProcess.smk --configfile ./config.yml --cores ${limited} --latency-wait 120
    cat 3readlength/*.txt > Length_count.txt
    ##
    rmp_extract.py --path ./4rpmData/ --total $mirsum --least $mirleast --sample $mirsamnum --output ./5mirdpOut/reads_18_26.fa --log rmp_log.txt
    cp ./5mirdpOut/reads_18_26.fa reads_18_26.fa
    ### read alignment information
    bowtie -p ${limited} -v 0 -f -t -a -m $mirmap --un ./5mirdpOut/unreads_18_26.fa \
        ${species_path}/Genome/Genome ./5mirdpOut/reads_18_26.fa >./5mirdpOut/mapping_location.sam
    bowtie -p ${limited} -v 1 -f -t -m $mirmap -a ${species_path}/Genome/Genome ./5mirdpOut/unreads_18_26.fa \
    >./5mirdpOut/mismapping.sam

    cat ./5mirdpOut/mapping_location.sam ./5mirdpOut/mismapping.sam | \
    awk '$1!~/@/{OFS="\t";if(NF==8){print $1,$4,$5+1,$5+length($6),$3,$2,"-"}else{print $1,$4,$5+1,$5+length($6),$3,$2,$9}}' - >./5mirdpOut/reads_18_26.txt
    rm ./5mirdpOut/mapping_location.sam ./5mirdpOut/mismapping.sam
    cp ./5mirdpOut/reads_18_26.txt reads_18_26.txt

    ## Identification tool
    if [ $method == "mirdp" ];then
        sed "s/ /_x/" ./5mirdpOut/reads_18_26.fa >./5mirdpOut/reads_18_26.flat
        ### miRDP2
        if [ -d ./5mirdpOut/reads_18_26-${mirmap}-0-${mirleast} ];then
            rm -r ./5mirdpOut/reads_18_26-${mirmap}-0-${mirleast}
        fi
        bash ${script_path}/miRDP2/miRDP2_pipeline.bash \
            -L ${mirmap} -M 0 -R ${mirleast} \
            -d ${script_path}/miRDP2 \
            -g ${species_path}/Genome/Genome.fa \
            -i ${species_path}/Genome/Genome \
            -f ./5mirdpOut/reads_18_26.flat -o ./5mirdpOut

        cp ./5mirdpOut/reads_18_26-${mirmap}-0-${mirleast}/*_P_prediction ./5mirdpOut/prediction.txt

        python ${script_path}/parsemiRDP2.py --fasta ./5mirdpOut/reads_18_26.fa \
            --genome ${species_path}/Genome/Genome.fa --faifile ${species_path}/Genome/Genome.fa.fai \
            --txtfile ./5mirdpOut/prediction.txt --output ./5mirdpOut/out_info.gff3

        python ${script_path}/parseMirGFF.py --input ./5mirdpOut/out_info.gff3 --fasta ./5mirdpOut/reads_18_26.fa --output  ./5mirdpOut/out_info.txt
        cp ./5mirdpOut/out_info.txt $outxt
    elif [ $method == "mircat" ];then
        awk -v RS=">" 'NR>1{print ">"$1"("$2")\n"$3}' 5mirdpOut/reads_18_26.fa >4mircatOut/mircat_reads_18_26.fa
        ### setting of mircat2
        echo """{
    "\"srna_files\"": [
        {
            "\"srna_filename\"": "\"${out_dir}/4mircatOut/mircat_reads_18_26.fa\""
        }
    ],
    "\"genome_filename\""    : "\"${species_path}/Genome/Genome.fa\"",
    "\"miRCAT2_Output_Dir\"" : "\"${out_dir}/4mircatOut\"",
    "\"miRCAT2_params\""     : "\"${script_path}/UEA_Workbench/default_miRCat2_plant_params.cfg\"",
    "\"annotation_filename\"": "\"${species_path}/GFF3/Annotation.gff3\""
}
""" >./default_mircat2.json
        ### mircat2
        java -jar -Xmx16g ${script_path}/UEA_Workbench/ServerWorkbench.jar \
        -tool mircat2 -config ./default_mircat2.json -standalone

        python ${script_path}/parsemircat.py --fasta ./4mircatOut/mircat_reads_18_26.fa \
            --genome ${species_path}/Genome/Genome.fa --faifile ${species_path}/Genome/Genome.fa.fai \
            --txtfile ./4mircatOut/mircat_reads_18_26_output.csv --output ./4mircatOut/out_info.txt
        cp ./4mircatOut/out_info.txt $outxt
    fi
    ### end identification
    ### output data
    cp ${curdir}/scripts/srnaProcess.Rmd ${out_dir}
    Rscript -e "rmarkdown::render('${out_dir}/srnaProcess.Rmd')"
    mkdir -p ${outputhtml%.*}_files
    mv srnaProcess_files ${outputhtml%.*}_files/srnaProcess_files
    cp srnaProcess.html ${outputhtml}
    cp 1multiQC/multiQC.html ${multiqcres}
    tar -zcvf $exptar -C ./4rpmData/ .
}

main "$@"
