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
        
rule map_reads:
    input:
        reads=get_trimmed_reads
    output:
        temp("%s/mapped/{sample}.sorted.bam" % (config["project-folder"]))
    log:
        "%s/logs/bwa_mem/{sample}.log" % (config["project-folder"])
    params:
        index=config["ref"]["genome"],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate"
    threads: 8
    wrapper:
        "0.27.1/bio/bwa/mem"


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
#        bai=get_recal_input(bai=True),
        ref=config["ref"]["genome"],
#        known=config["ref"]["known-variants"]
    output:
        bam=protected("%s/recal/{sample}.bam" % (config["project-folder"]))
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"]
    log:
        "%s/logs/gatk/bqsr/{sample}.log" % (config["project-folder"])
    wrapper:
        "0.27.1/bio/gatk/baserecalibrator"


rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    wrapper:
        "0.27.1/bio/samtools/index"
