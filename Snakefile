configfile: "config.yaml"

rule all:
    input:
        "TestData/" + "_".join(config["sra"]) + ".spritz.tr.snpeff.protein.withmods.xml"

include: "rules/downloads.smk"
include: "rules/align.smk"
include: "rules/variants.smk"
include: "rules/isoforms.smk"
include: "rules/proteogenomics.smk"
