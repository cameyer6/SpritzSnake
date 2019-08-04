FUSION_REF_VERSION = GENOME_VERSION + "_gencode_v29_CTAT_lib_Mar272019.plug-n-play"

rule download_premade_fusion_indices:
    '''Get the premade STAR-Fusion indices'''
    output: "data/" + FUSION_REF_VERSION + "/ctat_genome_lib_build_dir/ref_genome.fa"
    shell:
        "wget -O - https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/" + FUSION_REF_VERSION + ".tar.gz | "
        "tar -C data -xz"

rule rsem_star_fusion:
    '''Align with chimeric alignments and analyze coding effects with STAR-Fusion'''
    input:
        tmpdir=temp(directory("tmp")),
        genomelibsa="data/" + FUSION_REF_VERSION + "/ctat_genome_lib_build_dir/ref_genome.fa",
        genomelibdir=directory("data/" + FUSION_REF_VERSION + "/ctat_genome_lib_build_dir"),
        fq1="data/trimmed/{sra}.trim_1.fastq.gz" if check_sra() is True else expand("data/{fq1}_1.fastq.gz", fq1=config["fq1"]),
        fq2="data/trimmed/{sra}.trim_2.fastq.gz" if check_sra() is True else expand("data/{fq2}_2.fastq.gz", fq2=config["fq2"]),
    output:
        "output/{sra}FusionAnalysis/star-fusion.fusion_predictions.abridged.coding_effect.tsv"
    resources: mem_mb=50000
    threads: 12
    log: "output/{sra}STARFusion.log"
    shell:
        "(STAR-Fusion --examine_coding_effect --CPU {threads} --tmpdir {input.tmpdir} "
        " --genome_lib_dir {input.genomelibdir} --output_dir output/{sra}FusionAnalysis "
        " --left_fq <(zcat {input.fq1}) --right_fq <(zcat {input.fq2})) &> {log}"

rule generate_fusion_proteins:
    '''Use coding effects to generate fusion proteins'''
    input:
        expand("output/{sra}FusionAnalysis/star-fusion.fusion_predictions.abridged.coding_effect.tsv", sra=config["sra"]),
        transfermods=TRANSFER_MOD_DLL,
    output:
        "output/FusionProteins.xml",
        "output/FusionProteins.withmods.xml"
    shell:
        "dotnet {input.transfermods} -x {input.unixml} -f " +
        ",".join(expand("output/{sra}FusionAnalysis/star-fusion.fusion_predictions.abridged.coding_effect.tsv", sra=config["sra"]))
