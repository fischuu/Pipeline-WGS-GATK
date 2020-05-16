import pandas as pd
from snakemake.utils import validate
from multiprocessing import cpu_count

report: "report/workflow.rst"

###### Sample sheets #####
#samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
sampleList= pd.read_table(config["samples"], header=None)[0].tolist()

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
#units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
#validate(units, schema="schemas/units.schema.yaml")

# contigs in reference genome
contigs = pd.read_table(config["ref"]["genome"] + ".fai",
                        header=None, usecols=[0], squeeze=True, dtype=str)

##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(sampleList)

##### Helper functions #####

def get_fastq(wildcards):
    """Get fastq files of given sample."""
    fastqs = units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
    return {"r1": fastqs.fq1, "r2": fastqs.fq2}


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample), "platform"])


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample."""
    return expand("%s/trimmed/{sample}.{group}.fastq.gz" % (config["project-folder"]),
                      group=[1, 2], **wildcards)

def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("%s/recal/{sample}.bam" % (config["project-folder"]),
                  sample=wildcards.sample)


def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


def get_call_variants_params(wildcards, input):
    return (get_regions_param(regions=input.regions, default=f"--intervals {wildcards.contig}") +
            config["params"]["gatk"]["HaplotypeCaller"])


def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = "%s/mapped/{sample}.sorted.bam" % (config["project-folder"])
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "%s/dedup/{sample}.bam" % (config["project-folder"])
    if bai:
        if config["processing"].get("restrict-regions"):
            # case 3: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 4: no index needed
            return []
    else:
        return f


##### Target rules #####

rule all:
    input:
        "%s/annotated/all.vcf.gz" % (config["project-folder"]),
        "%s/qc/multiqc.html" % (config["project-folder"]),
        "%s/plots/depths.svg" % (config["project-folder"]),
        "%s/plots/allele-freqs.svg" % (config["project-folder"])


##### Modules #####

include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/stats.smk"
include: "rules/qc.smk"
include: "rules/annotation.smk"
