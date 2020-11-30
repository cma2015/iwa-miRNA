#!/bin/bash
# Genome_version_change.sh

set -x
script_name="verchange.sh"
script_ver="1.0.0"

# Help function
usage() {
  echo "--path     --file path"
  echo "--genome   --Genome file"
  echo "--flatrna  --flat format of miRNA annotation"
  echo "--version  --Version of script"
  echo "Example: bash ~/sRNAbox/tools/srnaDatabase/scripts/verchange.sh \
--path ~/sRNAbox/tools/srnaDatabase/tmp/20200112 \
--genome ~/sRNAbox/tools/srnaVis/Index/release-43/Zea_mays.B73_RefGen_v4.43/ \
--flatrna ~/sRNAbox/tools/srnaDatabase/tmp/20200112/PmiREN.txt"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){
    ## Parsing options
    ARGS=`getopt -o vh --long version,help,path:,genome:,flatrna: -n "$0" -- "$@"`
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
            --genome) genome=$2; shift 2;;
            --flatrna) flatrna=$2; shift 2;;
            -v|--version) version;;
            -h|--help) usage;;
            --) shift; break;;
            *)  echo "Internal error!"; exit 1;;
        esac
    done

    mkdir -p ${curdir}/GMAP
    mv ${flatrna} ${curdir}/GMAP/tmp.txt
    awk 'NR>1{print ">"$1"\n"$3}' ${curdir}/GMAP/tmp.txt >${curdir}/GMAP/tmp.fa

    if [ ! -d ${genome}/gmap_index ];then
        gmap_build -D ${genome} -d gmap_index -k 15 ${genome}/Genome/Genome.fa
    fi

    gmap -D ${genome} -d gmap_index --nosplicing --no-chimeras \
    --allow-close-indels=0 -n 10 -f psl -t 20 ${curdir}/GMAP/tmp.fa | awk '$2==0&&$18==1' > ${curdir}/GMAP/tmp.psl
    python ${curdir}/gmap_change.py --inputfile ${curdir}/GMAP/tmp.txt --psl ${curdir}/GMAP/tmp.psl --outfile ${flatrna}
}

main "$@"
