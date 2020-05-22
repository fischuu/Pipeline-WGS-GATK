rule trim_reads:
    input:
        unpack(get_fastq)
    output: 
        r1=temp("%s/trimmed/{sample}.1.fastq.gz") % (config["project-folder"]),
        r2=temp("%s/trimmed/{sample}.2.fastq.gz") % (config["project-folder"]),
        unpaired="%s/trimmed/{sample}.final.pe_se.fastq.gz"% (config["project-folder"]),
        trimlog="%s/trimmed/{sample}.trimlog.txt" % (config["project-folder"])
    params:
        r1_unpaired=temp("%s/trimmed/{sample}.1.unpaired.fastq.gz") % (config["project-folder"]),
        r2_unpaired=temp("%s/trimmed/{sample}.2.unpaired.fastq.gz") % (config["project-folder"])
    log:
        "%s/logs/trimmomatic/{sample}.log" % (config["project-folder"])
    threads:
        lambda cores: cpu_count()
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        trimmomatic PE \
            -threads {threads} \
            -trimlog {output.trimlog} \
            <(gzip -dc {input.r1} ) \
            <(gzip -dc {input.r2} ) \
            >(cut -f 1 -d " " | gzip -9 > {output.r1} ) \
            {params.r1_unpaired} \
            >(cut -f 1 -d " " | gzip -9 > {output.r2} ) \
            {params.r2_unpaired} \
            ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
            2> {log}
            
        zcat {params.r1_unpaired} {params.r2_unpaired} |
        cut -f 1 -d " " |
        gzip -9 > {output.unpaired}
        
        rm {params.r1_unpaired} {params.r2_unpaired}
        """

rule bwa_create_index:
        """
        Index the Reference Genome (BWA).
        """
        input:
            "%s" % (config["ref"]["genome"])
        output:
            "%s" % (config["ref"]["genome-bwa-index"])
        log:
            "%s/logs/bwa-IndexReferenceGenome.log" % (config["project-folder"])
        benchmark:
            "%s/benchmark/bwaIndexReferenceGenome.benchmark.tsv" % (config["project-folder"])
        threads: lambda cores: cpu_count()
        shell:"""
            echo "Number of threads used:" {threads}
            bwa index -a bwtsw {input} 2> {log}
                samtools faidx {input} 2> {log}
      	"""  
        
rule map_reads:
    """
    Map and sort the alignment (BWA).
    """
    input:
        reads=get_trimmed_reads,
        index="%s" % (config["ref"]["genome-bwa-index"])
    output:
        "%s/mapped/{sample}.sam" % (config["project-folder"])
    log:
        "%s/logs/bwa_mem/{sample}.log" % (config["project-folder"])
    params:
        index=config["ref"]["genome"],
        extra=get_read_group
    threads: lambda cores: cpu_count()
    conda: "../envs/mapping.yaml"
    shell:"""
        bwa mem -t {threads} {params.extra} {params.index} {input.reads} > {output} 2> {log}
    """
    
rule sort_alignment:
    """
    Map and sort the alignment (BWA).
    """
    input:
       "%s/mapped/{sample}.sam" % (config["project-folder"])
    output:
        "%s/mapped/{sample}.sorted.bam" % (config["project-folder"])
    log:
        "%s/logs/samtools/{sample}.log" % (config["project-folder"])
    threads: lambda cores: cpu_count()
    conda: "../envs/mapping.yaml"
    shell:"""
        samtools sort {input} > {output} 2> {log}
    """    

rule mark_duplicates:
    input:
        "%s/mapped/{sample}.sorted.bam" % (config["project-folder"])
    output:
        bam=temp("%s/dedup/{sample}.bam" % (config["project-folder"])) ,
        metrics="%s/qc/dedup/{sample}.metrics.txt" % (config["project-folder"])
    log:
        "%s/logs/picard/dedup/{sample}.log" % (config["project-folder"])
    params:
        config["params"]["picard"]["MarkDuplicates"]
    wrapper:
        "0.26.1/bio/picard/markduplicates"


rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
#        bai=[],#get_recal_input(bai=True),
#        ref=config["ref"]["genome"],
#        known=[]
    output:
        bam=protected("%s/recal/{sample}.bam" % (config["project-folder"]))
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"]
    log:
        "%s/logs/gatk/bqsr/{sample}.log" % (config["project-folder"])
#    wrapper:
#        "0.27.1/bio/gatk/baserecalibrator"
    shell:"""
      cp {input.bam} {output}
    """
    
rule index_recalibrated_bams:
    input:
        bam=get_sample_bams
    output:
        bai="%s/recal/{sample}.bam.bai" % (config["project-folder"])
    log:
        "%s/logs/gatk/bqsr/{sample}.log" % (config["project-folder"])
    conda: "../envs/mapping.yaml"
    shell:"""
      samtools index {input.bam} 2> {log}
    """    