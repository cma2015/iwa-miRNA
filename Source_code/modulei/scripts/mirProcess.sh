#!/bin/bash
# mirProcess.sh

set -x
script_name="mirProcess.sh"
script_ver="1.0.0"

# Help function
usage() {
  echo "--path     --file path"
  echo "--fastagz  --fasta path"
  echo "--species  --species path"
  echo "--filter   --limited sample names"
  echo "--readnum  --read number"
  echo "--min      --min length"
  echo "--max      --max length"
  echo "--map      --multile alignment"
  echo "--sum      --total RPM"
  echo "--least    --at least RPM"
  echo "--samnum   --at least samples"
  echo "--method   --method using for identification"
  echo "--thread   --Number of threads"
  echo "--outxt    --miRDP2 results"
  echo "--outfasta --sequences"
  echo "--outmap   --alignment file"
  echo "--expfile  --result path"   
  echo "--version  --Version of script"
  echo "Example: bash ~/sRNAbox/tools/modulei/scripts/mirProcess.sh --path ~/sRNAbox/tools/modulei \
--fastagz ~/sRNAbox/tools/modulei/tmp/fasta.tar.gz --species ~/sRNAbox/tools/modulei/Index/arabidopsis_thaliana_43 \
--readnum 1000000 --min 20 --max 24 --map 50 --sum 50 --least 5 --samnum 2 --method mirdp --thread 1 \
--outxt ~/out_info.txt --outfasta ~/reads_18_26.fa --outmap ~/reads_18_26.fa --expfile ~/path.out"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){
    ## Parsing options
    ARGS=`getopt -o vh --long version,help,path:,fastagz:,species:,:filter,readnum:,\
min:,max:,map:,sum:,least:,samnum:,method:,thread:,outxt:,outfasta:,outmap:,expfile: -n "$0" -- "$@"`
    if [ $? != 0 ]; then
        echo "Terminating..."
        exit 1
    fi
    eval set -- "${ARGS}"
    echo formatted parameters=[$@]
    while true
    do
        case "$1" in
            --path) curdir=$2; shift 2;;
            --fastagz) fastagz=$2; shift 2;;
            --species) species=$2; shift 2;;
            --filter) filter=$2; shift 2;;
            --readnum) readnum=$2; shift 2;;
            --min) min=$2; shift 2;;
            --max) max=$2; shift 2;;
            --map) multimap=$2; shift 2;;
            --sum) sum=$2; shift 2;;
            --least) least=$2; shift 2;;
            --samnum) samnum=$2; shift 2;;
            --method) method=$2; shift 2;;
            --thread) thread=$2; shift 2;;
            --outxt) outxt=$2; shift 2;;
            --outfasta) outfasta=$2; shift 2;;
            --outmap) outmap=$2; shift 2;;
            --expfile) expfile=$2; shift 2;;
            -v|--version) version;;
            -h|--help) usage;;
            --) shift; break;;
            *)  echo "Internal error!"; exit 1;;
        esac
    done

    ## Add environment path
    scriptpath=${curdir}/scripts
    PATH="$PATH:${scriptpath}"

    # Define the output directory
    OUTDIR=${curdir}/tmp/111111 #2019_`date +%s%N | cut -c6-13`

    if [ ! -d $OUTDIR ]; then
        mkdir $OUTDIR $OUTDIR/0collapseData
        tar -zxvf $fastagz -C $OUTDIR/0collapseData
    fi
    cd ${OUTDIR}

    speciesPath=${curdir}/`cat $species`
    ## Data processing
    echo -e "inputpath:\n    $curdir:" >./config.yml
    echo -e "speciespath:\n    $speciesPath:" >>./config.yml
    echo -e "readnum:\n    '$readnum:'" >>./config.yml
    echo -e "multimap:\n    '$multimap':" >>./config.yml
    echo -e "sum:\n    '$sum':" >>./config.yml
    echo -e "least:\n    '$least':" >>./config.yml
    echo -e "samnum:\n    '$samnum':" >>./config.yml
    echo -e "outpath:\n    $OUTDIR:\nsamples:" >>./config.yml

    if [ -z $filter ];then
        ls $OUTDIR/0collapseData | sed -r 's/.+\///;s/\..+//' | while read line
        do
            echo "    ${line}:" >>./config.yml
        done
    else
        cat $filter | while read line
        do
            echo "    ${line}:" >>./config.yml
        done
    fi

    if [ ! -d $OUTDIR/2rpmData ]; then    
    snakemake --snakefile ${scriptpath}/mirProcess.smk --configfile ./config.yml --cores ${thread} --latency-wait 120
    fi

    mkdir -p ./3mirdpOut
    rmp_extract.py --path ./2rpmData/ --total $sum --least $least --sample $samnum --output ./3mirdpOut/reads_18_26.fa --log rmp_log.txt

    ## Identification tool
    if [ $method == "mirdp" ];then
        
        sed "s/ /_x/" ./3mirdpOut/reads_18_26.fa >./3mirdpOut/reads_18_26.flat
        ### miRDP2
        if [ ! -d ./3mirdpOut/reads_18_26-20-0-5 ];then
            bash ${scriptpath}/miRDP2/miRDP2_pipeline.bash \
                -L 20 -M 0 -R 5 \
                -d ${scriptpath}/miRDP2 \
                -g ${speciesPath}/Genome/Genome.fa \
                -i ${speciesPath}/Genome/Genome \
                -f ./3mirdpOut/reads_18_26.flat -o ./3mirdpOut
        fi

        cp ./3mirdpOut/reads_18_26-20-0-5/*_P_prediction ./3mirdpOut/prediction.txt
        cp ./3mirdpOut/reads_18_26.fa ./3mirdpOut/reads_18_26.fa

        ### read alignment information
        bowtie -p 5 -v 0 -f -t -a -m $multimap --un ./3mirdpOut/unreads_18_26.fa \
            ${speciesPath}/Genome/Genome ./3mirdpOut/reads_18_26.fa >./3mirdpOut/mapping_location.sam
        bowtie -p 5 -v 1 -f -t -m $multimap -a ${speciesPath}/Genome/Genome ./3mirdpOut/unreads_18_26.fa \
        >./3mirdpOut/mismapping.sam
        
        cat ./3mirdpOut/mapping_location.sam ./3mirdpOut/mismapping.sam | \
        awk '$1!~/@/{OFS="\t";if(NF==8){print $1,$4,$5+1,$5+length($6),$3,$2,"-"}else{print $1,$4,$5+1,$5+length($6),$3,$2,$9}}' - >./3mirdpOut/reads_18_26.txt
        rm ./3mirdpOut/mapping_location.sam ./3mirdpOut/mismapping.sam
        
        python ${scriptpath}/parsemiRDP2.py --fasta ./3mirdpOut/reads_18_26.fa \
            --genome ${speciesPath}/Genome/Genome.fa --faifile ${speciesPath}/Genome/Genome.fa.fai \
            --txtfile ./3mirdpOut/prediction.txt --output ./out_info.gff3
        # RESULTDIR=./3mirdpOut

        python ${scriptpath}/parseMirGFF.py --input ./out_info.gff3 --output  ./out_info.txt
        
        ## GFF3 file and miRNAs plot
        cp ./out_info.txt $outxt
        cp ./3mirdpOut/reads_18_26.fa $outfasta
        cp ./3mirdpOut/reads_18_26.txt $outmap
        echo "${OUTDIR}" >$expfile
        tar -zcvf $expfile -C ${OUTDIR}/2rpmData .

    elif [ $method == "mircat" ];then
        ### setting of mircat2
        echo """
{
    "srna_files": [
        {
            "srna_filename": "./3mirdpOut/reads_18_26.fa"
        }
    ],
    "genome_filename"    : "${species}/Genome/Genome.fa",
    "miRCAT2_Output_Dir" : "./4mircatOut",
    "miRCAT2_params"     : "${scriptpath}/UEA_Workbench/default_miRCat2_plant_params.cfg",
    "annotation_filename": "${species}/GFF3/Annotation.gff3"
}
""" >>./default_mircat2.json
        ### mircat2
        java -jar -Xmx16g ${scriptpath}/UEA_Workbench/ServerWorkbench.jar \
        -tool mircat2 -config ./default_mircat2.json -standalone

    fi

}

main "$@"
