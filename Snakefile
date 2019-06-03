configfile: "config.yaml"

def check_dir():
    if ('analysisDirectory' in config and config["analysisDirectory"] is not None) and len(config["analysisDirectory"]) > 0:
            return True
    return False

def output(wildcards):
     if check_dir():
         return expand(["{dir}/combined.spritz.snpeff.protein.withmods.xml", "{dir}/combined.spritz.isoform.protein.withmods.xml", "{dir}/GRCh38.86.protein.withmods.xml", "{dir}/combined.spritz.isoformvariants.protein.withmods.xml"], dir=config["analysisDirectory"])
     return expand(["output/combined.spritz.snpeff.protein.withmods.xml", "output/combined.spritz.isoform.protein.withmods.xml", "output/GRCh38.86.protein.withmods.xml", "output/combined.spritz.isoformvariants.protein.withmods.xml"])

rule all:
    input: output

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
