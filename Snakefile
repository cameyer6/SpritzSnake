configfile: "config.yaml"

rule all:
    input:
        "data/combined.spritz.snpeff.protein.withmods.xml.gz",
        "data/combined.spritz.noindels.snpeff.protein.withmods.xml.gz",
        "data/combined.spritz.isoformvariants.protein.withmods.xml.gz",
        "data/combined.spritz.noindels.isoformvariants.protein.withmods.xml.gz",
        "data/combined.spritz.isoform.protein.withmods.xml.gz", # no variants
        "data/GRCh38.86.protein.withmods.xml.gz", # no variants
        # "clean_snpeff"

rule clean:
    shell:
        "rm -rf data/ ChromosomeMappings/ SnpEff/ tmp/ fast.tmp/ && "
        "cd GtfSharp && dotnet clean && cd .. && "
        "cd TransferUniProtModifications && dotnet clean && cd .."

include: "rules/downloads.smk"
include: "rules/align.smk"
include: "rules/variants.smk"
include: "rules/isoforms.smk"
include: "rules/proteogenomics.smk"
include: "rules/qc.smk"
