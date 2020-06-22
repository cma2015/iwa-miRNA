#! bin/bash

set -x
script_name="obtainSpecies.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "--curpath         --file path"
  echo "--ensembl         --species name in ensembl plant"
  echo "--verorder        --species verorder in ensembl plant"
  echo "--outgenome       --species verorder in ensembl plant"
  echo "--outgff3         --species verorder in ensembl plant"
  echo "--outtranscript   --species verorder in ensembl plant"
  echo "--outtrsnsnoRNA   --species verorder in ensembl plant"
  echo "--name            --species name"
  echo "--genome          --input genome"
  echo "--gff             --input gff"
  echo "--transcript      --input transcript"
  echo "--otherrna        --input otherrna"
  echo "--speciesout      --output species path"
  echo "--verorder        --verorder of script"
  echo "Example: ./$script_name "
  exit 1
}
# verorder function
verorder(){
    echo "$script_name $script_ver"
    exit 1
}
main(){
    # Parsing options
    ARGS=`getopt -o vh --long version,help,curpath:,ensembl:,verorder:,outgenome:,outgff3:,outtranscript:,outtrsnsnoRNA:,name:,genome:,gff:,transcript:,otherrna:,speciesout: -n "$0" -- "$@"`
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
            --ensembl) ensembl=$2; shift 2;;
            --verorder) verorder=$2; shift 2;;
            --outgenome) outgenome=$2; shift 2;;
            --outgff3) outgff3=$2; shift 2;;
            --outtranscript) outtranscript=$2; shift 2;;
            --outtrsnsnoRNA) outtrsnsnoRNA=$2; shift 2;;
            --name) speciesname=$2; shift 2;;
            --genome) genome=$2; shift 2;;
            --gff) gff=$2; shift 2;;
            --transcript) transcript=$2; shift 2;;
            --otherrna) otherrna=$2; shift 2;;
            --speciesout) speciesout=$2; shift 2;;
            -v|--version) version;;
            -h|--help) usage;;
            --) shift; break;;
            *)  echo "Internal error!"; exit 1;;
        esac
    done

    if [ -z $ensembl ];then
    	# Define the output directory
        out_dir=${curdir}/Index/${speciesname}_`date +%s%N | cut -c6-8`
        if [ ! -d $out_dir ]; then
            mkdir $out_dir $out_dir/Genome $out_dir/GFF3 \
            $out_dir/Transcripts $out_dir/trsnsnoRNAs
        fi
        if [ -f $genome ];then
            mv $genome ${out_dir}/Genome/Genome.fa && \
            bowtie-build ${out_dir}/Genome/Genome.fa ${out_dir}/Genome/Genome
        else
            wget $genome -O ${out_dir}/Genome/Genome.fa && \
            bowtie-build ${out_dir}/Genome/Genome.fa ${out_dir}/Genome/Genome
        fi
        if [ -f $otherrna ];then
            mv $otherrna ${out_dir}/trsnsnoRNAs/trsnsnoRNAs.fa && \
            bowtie-build ${out_dir}/trsnsnoRNAs/trsnsnoRNAs.fa ${out_dir}/trsnsnoRNAs/trsnsnoRNAs
        else
            wget $genome -O ${out_dir}/trsnsnoRNAs/trsnsnoRNAs.fa && \
            bowtie-build ${out_dir}/trsnsnoRNAs/trsnsnoRNAs.fa ${out_dir}/trsnsnoRNAs/trsnsnoRNAs
        fi
        if [ -f $gff ];then
            mv $gff ${out_dir}/GFF3/Annotation.gff3
        else
            wget $gff -O ${out_dir}/GFF3/Annotation.gff3
        fi
        if [ -f $transcript ];then
            mv $transcript $out_dir/Transcripts/transcripts.fa
        else
            wget $transcript -O $out_dir/Transcripts/transcripts.fa
        fi
        echo ${out_dir} >${speciesout}
    else
    	# Define the output directory
        out_dir=${curdir}/Index/${ensembl}_${verorder}
        if [ -d "${out_dir}" ];then
            echo "Index/${ensembl}_${verorder}" >${speciesout}
        else
            mkdir -p $out_dir
            cd $out_dir
            wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${verorder}/species_EnsemblPlants.txt -O EnsemblPlants.txt
            speciesFileName=`grep "${ensembl}" EnsemblPlants.txt | awk -F"\t" -v v=${verorder}  '{gsub(" ", "_",$5);print $2"."$5"."v}'`
            speciesGenome=`grep "${ensembl}" EnsemblPlants.txt | awk -F"\t" '{gsub(" ", "_",$5);print $2"."$5}'`
            #`grep "${ensembl}" EnsemblPlants.txt | awk -F"\t" '{gsub(" ", "_", $1);gsub(" ", "_",$5);print $1"."$5}'`
            mkdir Genome GFF3 Transcripts trsnsnoRNAs
            ## Genome file
            wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${verorder}/fasta/${ensembl}/dna/${speciesGenome^}.dna_sm.toplevel.fa.gz -O Genome/Genome.fa.gz && gzip -d Genome/Genome.fa.gz && bowtie-build Genome/Genome.fa Genome/Genome &
            ## GFF3 file
            wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${verorder}/gff3/${ensembl}/${speciesFileName^}.gff3.gz -O GFF3/Annotation.gff3.gz && gzip -d GFF3/Annotation.gff3.gz &
            ## Transcripts file
            wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${verorder}/fasta/${ensembl}/cdna/${speciesGenome^}.cdna.all.fa.gz -O Transcripts/transcripts.fa.gz && gzip -d Transcripts/transcripts.fa.gz &
            ## Proteins file
            # wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${verorder}/fasta/${ensembl}/pep/${speciesGenome}.pep.all.fa.gz -O Proteins/proteins.fa.gz && gzip -d Proteins/proteins.fa.gz &
            ## trsnsnoRNAs file
            wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${verorder}/fasta/${ensembl}/ncrna/${speciesGenome^}.ncrna.fa.gz -O trsnsnoRNAs/ncrna.all.fa.gz
            gzip -d trsnsnoRNAs/ncrna.all.fa.gz
            awk -v RS=">" '{if($0~/gene_biotype:snRNA|gene_biotype:snoRNA|gene_biotype:tRNA|gene_biotype:rRNA/){printf ">"$0}}' trsnsnoRNAs/ncrna.all.fa > trsnsnoRNAs/trsnsnoRNAs.fa
            bowtie-build trsnsnoRNAs/trsnsnoRNAs.fa trsnsnoRNAs/trsnsnoRNAs
            wait
            echo "Index/${ensembl}_${verorder}" >${speciesout}
        fi

        cd $out_dir
        
        if [ ! -z $outgenome ];then
            gzip -k Genome/Genome.fa && mv Genome/Genome.fa.gz $outgenome
        fi
        if [ ! -z $outgenome ];then
            cp GFF3/Annotation.gff3 $outgff3
        fi
        if [ ! -z $outgenome ];then
            cp Transcripts/transcripts.fa $outtranscript
        fi
        if [ ! -z $outgenome ];then
            cp trsnsnoRNAs/trsnsnoRNAs.fa $outtrsnsnoRNA
        fi
    fi
}
main "$@"
