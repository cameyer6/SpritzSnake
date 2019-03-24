rule all:
    input:
        "TestData/SRR393698.bam"

include: "rules/downloads.smk"
include: "rules/align.smk"
include: "rules/variants.smk"
