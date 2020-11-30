#! bin/bash
set -x
script_name="obtainSpecies.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "--curpath         --file path"
  echo "--ensembl         --species name in ensembl plant"
  echo "--verorder        --species verorder in ensembl plant"
  echo "--outgenome       --genome in ensembl plant"
  echo "--outgff3         --gff3 in ensembl plant"
  echo "--outtranscript   --transcript in ensembl plant"
  echo "--outprotein      --outprotein in ensembl plant"
  echo "--outtrsnsnoRNA   --species verorder in ensembl plant"
  echo "--name            --species name"
  echo "--genome          --input genome"
  echo "--gff             --input gff"
  echo "--transcript      --input transcript"
  echo "--protein         --input protein"
  echo "--otherrna        --input otherrna"
  echo "--speciesout      --output species path"
  echo "--version        --version of script"
  echo "Example: ./$script_name "
  exit 1
}

# version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){
    # Parsing options
    ARGS=`getopt -o vh --long curpath:,ensembl:,verorder:,outgenome:,outgff3:,outtranscript:,outtrsnsnoRNA:,outprotein:,infout:,\
name:,genome:,gff:,transcript:,protein:,otherrna:,speciesout:,randomname:,version,help -n 'example.bash' -- "$@"`
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
            --outprotein) outprotein=$2; shift 2;;
            --name) speciesname=$2; shift 2;;
            --genome) genome=$2; shift 2;;
            --gff) gff=$2; shift 2;;
            --transcript) transcript=$2; shift 2;;
            --protein) protein=$2; shift 2;;
            --otherrna) otherrna=$2; shift 2;;
            --speciesout) speciesout=$2; shift 2;;
            --randomname) randomname=$2; shift 2;;
            --infout) infout=$2; shift 2;;
            -v|--version) version;;
            -h|--help) usage;;
            --) shift; break;;
            *)  echo "Internal error!"; exit 1;;
        esac
    done

    if [ -z $ensembl ];then
        ##
        # Check for mandatory options
        # if [[ -f $genome ]] || [[ -f $gff ]] || [[ -f $otherrna ]] || [[ -f $transcript ]]; then
        #     echo "!!! Please check for inputs !!!"
        #     usage
        # fi

        # Define the output directory
        tmp_name=${speciesname/ /_}
        tmp_vec=`cat /dev/urandom | head -c 5 | md5sum | head -c 5`
        tmpfolder=${tmp_name}_${tmp_vec}
        out_dir=${curdir}/Index/${tmpfolder}

        if [ ! -d $out_dir ]; then
            mkdir $out_dir $out_dir/Genome $out_dir/GFF3 \
            $out_dir/Transcripts $out_dir/trsnsnoRNAs $out_dir/Proteins
        fi
        cd $out_dir
        # Upload or Download sequences and annotation file
        cp $genome Genome/Genome.fa
        bowtie-build Genome/Genome.fa Genome/Genome &
        cp $otherrna trsnsnoRNAs/trsnsnoRNAs.fa && \
        bowtie-build trsnsnoRNAs/trsnsnoRNAs.fa trsnsnoRNAs/trsnsnoRNAs
        cp $gff GFF3/Annotation.gff3 && cp $transcript Transcripts/transcripts.fa
        if [[ -f $protein ]]; then
            cp $protein Proteins/proteins.fa
        fi
        wait

        echo """
        ## The index construction have been completed.

        The path of ${speciesname} is: ${curdir}/Index/${tmpfolder}

        The directory contents were listed in a following tree structure:
        """ >${speciesout}
        tree -hfc | awk '{print "        "$0}' - >>${speciesout}
        python ${curdir}/scripts/genomeFeature.py ./
        Rscript ${curdir}/scripts/genomeInfo.R Genomic_info.pdf
        samtools faidx Genome/Genome.fa
        cp Genomic_info.pdf ${infout}
    else
        # Define the output directory
        tmp_name=${ensembl/ /_}
        if [ ${tmp_name}_${verorder} == "Arabidopsis_thaliana_47" ] || [ ${tmp_name}_${verorder} == "Zea_mays_47" ] || [ ${tmp_name}_${verorder} == "Triticum_aestivum_47" ] || \
        [ ${tmp_name}_${verorder} == "Oryza_sativa_47" ] || [ ${tmp_name}_${verorder} == "Glycine_max_47" ] || [ ${tmp_name}_${verorder} == "Gossypium_raimondii_47" ] ;then
            out_dir=${curdir}/Index/${tmp_name}_${verorder}
        else
            tmp_vec=`cat /dev/urandom | head -c 5 | md5sum | head -c 5`
            out_dir=${curdir}/Index/${tmp_name}_${verorder}_${tmp_vec}
        fi

        if [ ! -d "${out_dir}" ];then
            mkdir -p $out_dir
            cd $out_dir
            wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${verorder}/species_EnsemblPlants.txt -O EnsemblPlants.txt
            speciesFileName=`grep "${tmp_name,}" EnsemblPlants.txt | awk -F"\t" -v v=${verorder}  '{gsub(" ", "_",$5);print $2"."$5"."v}'`
            speciesGenome=`grep "${tmp_name,}" EnsemblPlants.txt | awk -F"\t" '{gsub(" ", "_",$5);print $2"."$5}'`
            #`grep "${ensembl}" EnsemblPlants.txt | awk -F"\t" '{gsub(" ", "_", $1);gsub(" ", "_",$5);print $1"."$5}'`
            mkdir Genome GFF3 Transcripts trsnsnoRNAs
            ## Genome file
            wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${verorder}/fasta/${tmp_name,}/dna/${speciesGenome^}.dna_sm.toplevel.fa.gz -O Genome/Genome.fa.gz
            gunzip Genome/Genome.fa.gz
            bowtie-build Genome/Genome.fa Genome/Genome &
            ## GFF3 file
            wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${verorder}/gff3/${tmp_name,}/${speciesFileName^}.gff3.gz -O GFF3/Annotation.gff3.gz
            gunzip GFF3/Annotation.gff3.gz
            ## Transcripts file
            wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${verorder}/fasta/${tmp_name,}/cdna/${speciesGenome^}.cdna.all.fa.gz -O Transcripts/transcripts.fa.gz
            gunzip Transcripts/transcripts.fa.gz
            # Proteins file
            wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${verorder}/fasta/${ensembl}/pep/${speciesGenome}.pep.all.fa.gz -O Proteins/proteins.fa.gz
            gunzip Proteins/proteins.fa.gz
            ## trsnsnoRNAs file
            wget ftp://ftp.ensemblgenomes.org/pub/plants/release-${verorder}/fasta/${tmp_name,}/ncrna/${speciesGenome^}.ncrna.fa.gz -O trsnsnoRNAs/ncrna.all.fa.gz
            gzip -d trsnsnoRNAs/ncrna.all.fa.gz
            awk -v RS=">" '{if($0~/gene_biotype:snRNA|gene_biotype:snoRNA|gene_biotype:tRNA|gene_biotype:rRNA/){printf ">"$0}}' trsnsnoRNAs/ncrna.all.fa > trsnsnoRNAs/trsnsnoRNAs.fa
            bowtie-build trsnsnoRNAs/trsnsnoRNAs.fa trsnsnoRNAs/trsnsnoRNAs
            wait
            samtools faidx Genome/Genome.fa
        fi

        cd $out_dir

        echo """
        ## The download of the reference sequence and the index construction have been completed.

        The path of ${ensembl} is: ${curdir}/Index/${tmp_name}_${verorder}

        The directory contents were listed in a following tree structure:
        """ >${speciesout}
        tree -hfc | awk '{print "        "$0}' - >>${speciesout}

        python ${curdir}/scripts/genomeFeature.py ./
        Rscript ${curdir}/scripts/genomeInfo.R Genomic_info.pdf

        if [ ! -f Genome/Genome.fa.gz ];then
            gzip -k Genome/Genome.fa && cp Genome/Genome.fa.gz $outgenome
        else
            cp Genome/Genome.fa.gz $outgenome
        fi
        cp GFF3/Annotation.gff3 $outgff3
        cp Transcripts/transcripts.fa $outtranscript
        cp trsnsnoRNAs/trsnsnoRNAs.fa $outtrsnsnoRNA
        if [ -f Proteins/proteins.fa ];then
            cp Proteins/proteins.fa $outprotein
        fi
        cp Genomic_info.pdf ${infout}
        # tmp_vec=`cat /dev/urandom | head -c 5 | md5sum | head -c 5`
        # echorandomname
    fi
}

main "$@"
