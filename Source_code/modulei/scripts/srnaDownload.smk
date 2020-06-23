
# Time: 2019-12
# sRNA-seq pipeline with snakemake

# configfile: "config.yml"

rule down_result:
    input:
        expand("{fileout}/{sample}.fastq", fileout=config["outpath"], sample=config["samples"])


rule download_sra:
    output: "{fileout}/{sample}.fastq"
    params:
        index="ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR",
        ipath=expand(config["inputpath"])
    shell:
        """
        set -x
        srrnam={wildcards.sample}
        if [[ -n $(echo $srrnam | grep "+") ]]; then
            srreach=${{srrnam//+/ }}
            for eachsam in $srreach
            do
                downurl=`curl https://www.ebi.ac.uk/ena/data/warehouse/filereport\?accession\=${{eachsam}}\&result\=read_run\&fields\=run_accession,fastq_ftp,fastq_md5,fastq_bytes | grep -oP "ftp.+?gz"`
                axel -n 10 $downurl
                gunzip $eachsam.fastq.gz
                cat *.fastq >>{output}
                rm *.fastq
            done
            gzip -k {output}
        else
            downurl=`curl https://www.ebi.ac.uk/ena/data/warehouse/filereport\?accession\=${{srrnam}}\&result\=read_run\&fields\=run_accession,fastq_ftp,fastq_md5,fastq_bytes | grep -oP "ftp.+?gz"`
            axel -n 10 $downurl
            gunzip $srrnam.fastq.gz
            mv $srrnam.fastq {output}
        fi
        gzip -k {output}
        """

# enaDataGet -f fastq $srrnam
# if [ -f $srrnam.fastq.gz ];then
#     gunzip $srrnam.fastq.gz
#     cat $srrnam.fastq >>{output}
# else
#     axel {params.index}/${{srrnam:0:6}}/${{srrnam}}/${{srrnam}}.sra
#     fastq-dump --split-3 -O {params.ipath} {params.ipath}/$srrnam.sra
#     cat {params.ipath}/$srrnam.fastq >>{output}
# fi
