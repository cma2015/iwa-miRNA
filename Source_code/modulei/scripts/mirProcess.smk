
# Time: 2019-12
# sRNA-seq pipeline with snakemake

rule sRNA_result:
    input:
        expand("{fileout}/2rpmData/{sample}.txt", fileout=config["outpath"], sample=config["samples"]),
        expand("{fileout}/2countData/{sample}.txt", fileout=config["outpath"], sample=config["samples"])
        
        

rule read_map:
    input: "{fileout}/0collapseData/{sample}.fasta"
    output:
        gm = "{fileout}/1mapData/{sample}_genome.fasta",
        clm = "{fileout}/1mapData/{sample}_rmtrsno.fasta"
    params:
        multimap=expand(config["multimap"]),
        speciespath=expand(config["speciespath"]),
        readnum=expand(config["readnum"])
    threads: 1
    shell:
        """
        set -x
        bowtie -p 5 -v 1 -f -t -m {params.multimap} -S --al {output.gm} {params.speciespath}/Genome/Genome {input} 1>/dev/null 2>&1
        bowtie -p 5 -v 0 -f -t -S --norc --un {output.clm} {params.speciespath}/trsnsnoRNAs/trsnsnoRNAs {output.gm} 1>/dev/null 2>&1
        clean=`awk '$1~/>/{{split($1,a,"-");b+=a[2]}}END{{print b}}' {input}`
        gg=`awk '$1~/>/{{split($1,a,"-");b+=a[2]}}END{{print b}}' {output.gm}`
        rmtr=`awk '$1~/>/{{split($1,a,"-");b+=a[2]}}END{{print b}}' {output.clm}`
        if [ $gg lt {params.readnum} ];then
            rm {wildcards.fileout}/1mapData/{wildcards.sample}_rmtrsno.fasta
        fi
        echo -e "{wildcards.sample}\t${{clean}}\t${{gg}}\t${{rmtr}}" >>{wildcards.fileout}/00Table_Summary_of_sRNA-seq.txt
        """

rule read_count:
    input: "{fileout}/1mapData/{sample}_rmtrsno.fasta"
    output: "{fileout}/2countData/{sample}.txt"
    threads: 1
    shell:
        """
        readmap=`awk '$1~/>/{{split($1,a,"-");b+=a[2]}}END{{print b}}' {input}`
        awk -v RS=">" -v SRR={wildcards.sample} 'NR>1{{split($1, a,"-");if(a[2]>5){{print SRR"\t"$2"\t"a[2]}}}}' {input} >{output}
        """

rule read_rpm:
    input: "{fileout}/1mapData/{sample}_rmtrsno.fasta"
    output: "{fileout}/2rpmData/{sample}.txt"
    threads: 1
    shell:
        """
        readmap=`awk '$1~/>/{{split($1,a,"-");b+=a[2]}}END{{print b}}' {input}`
        awk -v RS=">" -v SRR={wildcards.sample} 'NR>1{{split($1, a,"-");RPM=a[2]/'$readmap'*1000000;if(a[2]>5){{print SRR"\t"$2"\t"RPM}}}}' {input} >{output}
        """
