rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html="%s/qc/fastqc/{sample}.html" % (config["project-folder"]),
        zip="%s/qc/fastqc/{sample}.zip" % (config["project-folder"])
    wrapper:
        "0.27.1/bio/fastqc"


rule samtools_stats:
    input:
        "%s/recal/{sample}.bam" % (config["project-folder"])
    output:
        "%s/qc/samtools-stats/{sample}.txt" % (config["project-folder"])
    log:
        "%s/logs/samtools-stats/{sample}.log" % (config["project-folder"])
    wrapper:
        "0.27.1/bio/samtools/stats"

rule multiqc:
    input:
        expand(["%s/qc/samtools-stats/{sample}.txt" % (config["project-folder"]),
                "%s/qc/fastqc/{sample}.zip" % (config["project-folder"]),
                "%s/qc/dedup/{sample}.metrics.txt" % (config["project-folder"])],
               sample=sampleList),
        "%s/snpeff/all.csv" % (config["project-folder"])
    output:
        report("%s/qc/multiqc.html" % (config["project-folder"]), caption="../report/multiqc.rst", category="Quality control")
    log:
        "%s/logs/multiqc.log" % (config["project-folder"])
    wrapper:
        "0.27.1/bio/multiqc"
