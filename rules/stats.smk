rule vcf_to_tsv:
    input:
        "%s/annotated/all.vcf.gz" % (config["project-folder"])
    output:
        report("%s/tables/calls.tsv.gz" % (config["project-folder"]), caption="../report/calls.rst", category="Calls")
    conda:
        "../envs/rbt.yaml"
    shell:
        "bcftools view --apply-filters PASS --output-type u {input} | "
        "rbt vcf-to-txt -g --fmt DP AD --info ANN | "
        "gzip > {output}"


rule plot_stats:
    input:
        "%s/tables/calls.tsv.gz" % (config["project-folder"])
    output:
        depths=report("%s/plots/depths.svg" % (config["project-folder"]), caption="../report/depths.rst", category="Plots"),
        freqs=report("%s/plots/allele-freqs.svg" % (config["project-folder"]), caption="../report/freqs.rst", category="Plots")
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/plot-depths.py"
