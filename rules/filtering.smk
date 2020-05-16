def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL")


rule select_calls:
    input:
        ref=config["ref"]["genome"],
        vcf="%s/genotyped/all.vcf.gz" % (config["project-folder"])
    output:
        vcf=temp("%s/filtered/all.{vartype}.vcf.gz" % (config["project-folder"]))
    params:
        extra=get_vartype_arg
    log:
        "%s/logs/gatk/selectvariants/{vartype}.log" % (config["project-folder"])
    wrapper:
        "0.27.1/bio/gatk/selectvariants"


def get_filter(wildcards):
    return {
        "snv-hard-filter":
        config["filtering"]["hard"][wildcards.vartype]}


rule hard_filter_calls:
    input:
        ref=config["ref"]["genome"],
        vcf="%s/filtered/all.{vartype}.vcf.gz" % (config["project-folder"])
    output:
        vcf=temp("%s/filtered/all.{vartype}.hardfiltered.vcf.gz" % (config["project-folder"]))
    params:
        filters=get_filter
    log:
        "%s/logs/gatk/variantfiltration/{vartype}.log" % (config["project-folder"])
    wrapper:
        "0.27.1/bio/gatk/variantfiltration"


rule recalibrate_calls:
    input:
        vcf="%s/filtered/all.{vartype}.vcf.gz" % (config["project-folder"])
    output:
        vcf=temp("%s/filtered/all.{vartype}.recalibrated.vcf.gz" % (config["project-folder"]))
    params:
        extra=config["params"]["gatk"]["VariantRecalibrator"]
    log:
        "%s/logs/gatk/variantrecalibrator/{vartype}.log" % (config["project-folder"])
    wrapper:
        "0.27.1/bio/gatk/variantrecalibrator"


rule merge_calls:
    input:
        vcf=expand("%s/filtered/all.{vartype}.{filtertype}.vcf.gz" % (config["project-folder"]),
                   vartype=["snvs", "indels"],
                   filtertype="recalibrated"
                              if config["filtering"]["vqsr"]
                              else "hardfiltered")
    output:
        vcf="%s/filtered/all.vcf.gz" % (config["project-folder"])
    log:
        "%s/logs/picard/merge-filtered.log" % (config["project-folder"])
    wrapper:
        "0.27.1/bio/picard/mergevcfs"
