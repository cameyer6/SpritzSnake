rule all:
    input:
        "TestData/ERR315327_1.spritz.tr.snpeff.protein.withmods.xml"

include: "rules/downloads.smk"
include: "rules/align.smk"
include: "rules/variants.smk"
include: "rules/isoforms.smk"
include: "rules/proteogenomics.smk"
