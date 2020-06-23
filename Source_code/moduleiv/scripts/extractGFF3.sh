#!/bin/bash
# extractGFF3.sh

script_name="extractGFF3.sh"
script_ver="1.0.0"

# Help function
usage() {
  echo "--curpath    --file path"
  echo "--genegff    --Annotation of genome"
  echo "--TEanno     --TEanno"
  echo "--Geneanno   --Geneanno"
  echo "--outfile    --outfile"
  echo "--version    --Version of script"
  echo "Example: bash ~/sRNAbox/tools/moduleiv/scripts/extractGFF3.sh --curpath ~/sRNAbox/tools/moduleiv \
--genegff ~/sRNAbox/tools/moduleiv/tmp/Annotation.gff3 --TEanno ~/sRNAbox/tools/moduleiv/tmp/B73.structuralTEv2.filteredTE_Springer.gff3 \
--outfile out_test.tsv"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){
    ## Parsing options
    ARGS=`getopt -o vh --long version,help,curpath:,genegff:,TEanno:,Geneanno:,outfile: -n "$0" -- "$@"`
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
            --genegff) genegff=$2; shift 2;;
            --TEanno) TEanno=$2; shift 2;;
            --Geneanno) Geneanno=$2; shift 2;;
            --outfile) outfile=$2; shift 2;;
            -v|--version) version;;
            -h|--help) usage;;
            --) shift; break;;
            *)  echo "Internal error!"; exit 1;;
        esac
    done

    ## Add environment path
    SCRIPTPATH=${curdir}/scripts
    PATH="$PATH:${SCRIPTPATH}"

    # Define the output directory
    OUTDIR=${curdir}/tmp/20200110 #2019_`date +%s%N | cut -c6-13`
    if [ ! -d $OUTDIR ]; then
        mkdir -p ${curdir}/tmp $OUTDIR
    fi

    python ${SCRIPTPATH}/genomeFeature.py $genegff $OUTDIR/0gene.txt GENE

    if [[ ! -z $Geneanno ]]; then
        awk 'NR==FNR{a[$1]=$2;next}{OFS="\t";if(a[$4]){$NF=a[$4]};print $0}' $Geneanno $OUTDIR/0gene.txt > $OUTDIR/0result.txt
    else
        cp $OUTDIR/0gene.txt $OUTDIR/0result.txt
    fi

    if [[ ! -z $TEanno ]]; then
        python ${SCRIPTPATH}/genomeFeature.py $TEanno $OUTDIR/0TE.txt TE
        cat $OUTDIR/0result.txt $OUTDIR/0TE.txt | sort -k1,1 -k2,2n > ${outfile}
    else
        sort -k1,1 -k2,2n $OUTDIR/0result.txt > ${outfile}
    fi
}

main $@
