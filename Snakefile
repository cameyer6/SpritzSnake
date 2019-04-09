configfile: "config.yaml"

rule all:
    input:
        expand("TestData/{sample}_1.spritz.tr.snpeff.protein.withmods.xml", sample=config["sra"])

include: "rules/downloads.smk"
include: "rules/align.smk"
include: "rules/variants.smk"
include: "rules/isoforms.smk"
include: "rules/proteogenomics.smk"
