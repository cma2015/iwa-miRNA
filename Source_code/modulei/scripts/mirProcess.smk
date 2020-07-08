
# Time: 2019-12
# sRNA-seq pipeline with snakemake

rule sRNA_result:
    input:
        expand("{fileout}/3readlength/{sample}_all.txt", fileout=config["outpath"], sample=config["samples"]),
        expand("{fileout}/4rpmData/{sample}.txt", fileout=config["outpath"], sample=config["samples"]),
        expand("{fileout}/4countData/{sample}.txt", fileout=config["outpath"], sample=config["samples"])

rule read_map:
    input:
        fq =  "{fileout}/00rawdata/{sample}.fastq",
        cfq = "{fileout}/1cleandata/{sample}_clean.fastq",
        clfa = "{fileout}/1cleandata/{sample}_clean.fasta",
        cufa = "{fileout}/2collapsedata/{sample}.fasta"
    output:
        gm = "{fileout}/3mapData/{sample}_genome.fasta",
        clm = "{fileout}/3mapData/{sample}_rmtrsno.fasta",
        alltxt = "{fileout}/3readlength/{sample}_all.txt",
        unique = "{fileout}/3readlength/{sample}_unique.txt"
    params:
        multimap=expand(config["multimap"]),
        speciespath=expand(config["speciespath"]),
        readnum=expand(config["readnum"])
    threads: 1
    shell:
        """
        set -x
        bowtie -p 5 -v 1 -f -t -m {params.multimap} -S --al {output.gm} {params.speciespath}/Genome/Genome {input.cufa} 1>/dev/null 2>&1
        bowtie -p 5 -v 0 -f -t -S --norc --un {output.clm} {params.speciespath}/trsnsnoRNAs/trsnsnoRNAs {output.gm} 1>/dev/null 2>&1
        readnum=$((`cat {input.fq} | wc -l`/4))
        clean=$((`cat {input.cfq} | wc -l`/4))
        ratioval=`awk 'BEGIN{{printf ('$clean'/'$readnum')*100}}'`
        cleanlen=`grep ">" {input.clfa} | wc -l`
        gg=`awk '$1~/>/{{split($1,a,"-");b+=a[2]}}END{{print b}}' {output.gm}`
        rmtr=`awk '$1~/>/{{split($1,a,"-");b+=a[2]}}END{{print b}}' {output.clm}`
        rmtrunique=`grep ">" {output.clm} | wc -l`
        abundseq=`sed -n '2p' {output.clm}`
        abundnum=`sed -n '1p' {output.clm} | cut -d"-" -f 2`
        if [ $gg -lt {params.readnum} ];then
            rm {wildcards.fileout}/3mapData/{wildcards.sample}_rmtrsno.fasta
        fi
        echo -e "{wildcards.sample}\t${{readnum}}\t${{clean}}\t${{ratioval}}\t${{cleanlen}}\t${{gg}}\t${{rmtr}}\t${{rmtrunique}}\t${{abundseq}}\t${{abundnum}}" >>{wildcards.fileout}/00Table_Summary_of_sRNA-seq_data.txt
        awk -v sample={wildcards.sample} '$1!~/>/{{a[length($1)]++}}END{{OFS="\t";for(i in a)print "All",sample,i,a[i]}}' {input.clfa} >{output.alltxt}
        awk -v sample={wildcards.sample} '$1!~/>/{{a[length($1)]++}}END{{OFS="\t";for(i in a)print "Unique",sample,i,a[i]}}' {input.cufa} >{output.unique}
        """

rule read_count:
    input: "{fileout}/3mapData/{sample}_rmtrsno.fasta"
    output: "{fileout}/4countData/{sample}.txt"
    threads: 1
    shell:
        """
        readmap=`awk '$1~/>/{{split($1,a,"-");b+=a[2]}}END{{print b}}' {input}`
        awk -v RS=">" -v SRR={wildcards.sample} 'NR>1{{split($1, a,"-");if(a[2]>5){{print SRR"\t"$2"\t"a[2]}}}}' {input} >{output}
        """

rule read_rpm:
    input: "{fileout}/3mapData/{sample}_rmtrsno.fasta"
    output: "{fileout}/4rpmData/{sample}.txt"
    threads: 1
    shell:
        """
        readmap=`awk '$1~/>/{{split($1,a,"-");b+=a[2]}}END{{print b}}' {input}`
        awk -v RS=">" -v SRR={wildcards.sample} 'NR>1{{split($1, a,"-");RPM=a[2]/'$readmap'*1000000;if(a[2]>5){{print SRR"\t"$2"\t"RPM}}}}' {input} >{output}
        """
