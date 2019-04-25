configfile: "config.yaml"

rule all:
    input:
        ["data/combined.spritz.snpeff.protein.xml", "data/combined.spritz.isoformed.snpeff.protein.xml"]

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
