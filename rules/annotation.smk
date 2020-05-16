rule snpeff:
    input:
        "%s/filtered/all.vcf.gz" % (config["project-folder"]) 
    output:
        vcf=report("%s/annotated/all.vcf.gz" % (config["project-folder"]), caption="../report/vcf.rst", category="Calls"),
        csvstats="%s/snpeff/all.csv" % (config["project-folder"])
    log:
        "%s/logs/snpeff.log" % (config["project-folder"])
    params:
        reference=config["ref"]["name"],
        extra="-Xmx6g"
    wrapper:
        "0.27.1/bio/snpeff"
