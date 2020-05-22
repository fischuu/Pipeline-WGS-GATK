if "restrict-regions" in config["processing"]:
    rule compose_regions:
        input:
            config["processing"]["restrict-regions"]
        output:
            "%s/called/{contig}.regions.bed" % (config["project-folder"])
        conda:
            "../envs/bedops.yaml"
        shell:  
            "bedextract {wildcards.contig} {input} > {output}"


rule create_dict_file:
    input:
        config["ref"]["genome"]
    output:
        config["ref"]["genomeDict"]
    log:
        "%s/logs/gatk/createDict/createDict.log" % (config["project-folder"])
    conda: "../envs/calling.yaml"
    shell:"""
        java -jar CreateSequenceDictionary.jar R= {input} O= {output}
        """


rule call_variants:
    input:
        bam=get_sample_bams,
        ref=config["ref"]["genome"],
        dict=config["ref"]["genomeDict"],
        regions=[]
    output:
        gvcf=protected("%s/called/{sample}.g.vcf.gz" % (config["project-folder"]))
    log:
        "%s/logs/gatk/haplotypecaller/{sample}.log" % (config["project-folder"])
    params:
#        extra=get_call_variants_params
        extra=[]        
    wrapper:
        "0.27.1/bio/gatk/haplotypecaller"


rule combine_calls:
    input:
        ref=config["ref"]["genome"],
        gvcfs=expand("%s/called/{sample}.g.vcf.gz" % (config["project-folder"]), sample=sampleList)
    output:
        gvcf=temp("%s/called/all.g.vcf.gz" % (config["project-folder"]))
    log:
        "%s/logs/gatk/combinegvcfs.log" % (config["project-folder"])
    wrapper:
        "0.27.1/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        gvcf="%s/called/all.g.vcf.gz" % (config["project-folder"])
    output:
        vcf="%s/genotyped/all.vcf.gz" % (config["project-folder"])
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        "%s/logs/gatk/genotypegvcfs.log" % (config["project-folder"])
    wrapper:
        "0.27.1/bio/gatk/genotypegvcfs"


#rule merge_variants:
#    input:
#        vcf="%s/genotyped/all.vcf.gz" % (config["project-folder"])
#    output:
#        vcf="%s/genotyped/all_merged.vcf.gz" % (config["project-folder"])
#    log:
#        "%s/logs/picard/merge-genotyped.log" % (config["project-folder"])
#    wrapper:
#        "0.27.1/bio/picard/mergevcfs"
