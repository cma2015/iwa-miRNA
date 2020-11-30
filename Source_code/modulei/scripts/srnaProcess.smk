# Time: 2019-12
# sRNA-seq pipeline with snakemake

rule sRNA_result:
    input:
        #expand("{fileout}/1fastqc/{sample}_fastqc.html", fileout=config["outpath"], sample=config["samples"]),
        expand("{fileout}/2collapsedata/{sample}.fasta", fileout=config["outpath"], sample=config["samples"])

rule read_clean:
    input: "{fileout}/00rawdata/{sample}.fastq",
    output:
        cfastq = "{fileout}/1cleandata/{sample}_clean.fastq",
        ffasta = "{fileout}/1cleandata/{sample}.fasta"
    params:
        quality=expand(config["quality"]),
        minval=expand(config["min"]),
        maxval=expand(config["max"])
    threads: 1
    shell:
        """
        set -x
        tmpname={wildcards.fileout}/1cleandata/{wildcards.sample}
        adapter=`grep {wildcards.sample} {wildcards.fileout}/Adapter.txt |head -n 1| cut -f 2`
        phred=`grep {wildcards.sample} {wildcards.fileout}/Adapter.txt |head -n 1| cut -f 3`
        barcode=`grep {wildcards.sample} {wildcards.fileout}/Adapter.txt |head -n 1| cut -f 4`
        if [ $barcode == "-" ];then
            fastx_clipper -n -v -Q${{phred}} -l 10 -a ${{adapter}} -i {input} -o $tmpname.tmp1
        elif [ $barcode == "lf4" ];then
            fastx_clipper -n -v -Q${{phred}} -l 10 -a ${{adapter}} -i {input} | \
            fastx_trimmer -f 5 -i - | fastx_trimmer -t 4 -i - -o $tmpname.tmp1
        else
            fastx_clipper -n -v -Q${{phred}} -l 10 -a ${{adapter}} -i {input} | \
            fastx_trimmer -f $barcode -i - -o $tmpname.tmp1
        fi
        fastq_quality_filter -q {params.quality} -p 80 -Q${{phred}} -i $tmpname.tmp1 -o  {output.cfastq}
        rm $tmpname.tmp1
        fastq_to_fasta -r -n -v -Q${{phred}} -i {output.cfastq} -o {output.ffasta}
 
        mkdir -p {wildcards.fileout}/1fastqc
        fastqc -q -t 5 -O {wildcards.fileout}/1fastqc {input}
        fastqc -q -t 5 -O {wildcards.fileout}/1fastqc {output.cfastq}

        wc -l {input} | cut -d " " -f 1 >{wildcards.fileout}/00rawdata/{wildcards.sample}.fastq.count
        rm {input}
        wc -l {output.cfastq} | cut -d " " -f 1 >{wildcards.fileout}/1cleandata/{wildcards.sample}_clean.fastq.count
        """

rule length_filter:
    input: "{fileout}/1cleandata/{sample}.fasta"
    output: "{fileout}/1cleandata/{sample}_clean.fasta"
    threads: 1
    params:
        minval=expand(config["min"]),
        maxval=expand(config["max"])
    shell:
        """
        rm {wildcards.fileout}/1cleandata/{wildcards.sample}_clean.fastq
        sh /home/galaxy/tools/1_MiRNA_Compilation/scripts/length_cutoff.sh {input} {params.minval} {params.maxval} {output}
        """

rule collapse_fa:
    input: "{fileout}/1cleandata/{sample}_clean.fasta"
    output: "{fileout}/2collapsedata/{sample}.fasta"
    threads: 1
    shell:
        """
        fastx_collapser -v -i {input} -o {output}
        """

#rule fastqc:
#    input:
#        rawdata = "{fileout}/00rawdata/{sample}.fastq",
#        cleandata = "{fileout}/1cleandata/{sample}_clean.fastq"
#    output: "{fileout}/1fastqc/{sample}_fastqc.html"
#    threads: 1
#    shell:
#        """
#        set -x
#        mkdir -p {wildcards.fileout}/1fastqc
#        fastqc -q -t 5 -O {wildcards.fileout}/1fastqc {input.rawdata}
#        fastqc -q -t 5 -O {wildcards.fileout}/1fastqc {input.cleandata}
#        """
