rule all:
    input:
        "TestData/ERR315327_1.spritz.snpeff.vcf"

include: "rules/downloads.smk"
include: "rules/align.smk"
include: "rules/variants.smk"
include: "rules/isoforms.smk"
include: "rules/proteogenomics.smk"
