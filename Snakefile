rule all:
    input:
        "TestData/mapper0.bam"

include: "rules/downloads.smk"
include: "rules/align.smk"
include: "rules/variants.smk"
