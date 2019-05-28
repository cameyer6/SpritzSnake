UNIPROTXML="data/uniprot/human.protein.xml.gz" #"data/Homo_sapiens_202022.xml.gz"

rule download_protein_xml:
    output: UNIPROTXML
    shell: "python scripts/download_protein_xml.py | gzip -c > {output}"

rule build_transfer_mods:
    output: "TransferUniProtModifications/TransferUniProtModifications/bin/Release/netcoreapp2.1/TransferUniProtModifications.dll"
    log: "data/TransferUniProtModifications.build.log"
    shell:
        "(cd TransferUniProtModifications && "
        "dotnet restore && "
        "dotnet build -c Release TransferUniProtModifications.sln) &> {log}"

rule transfer_modifications1:
    input:
        transfermods="TransferUniProtModifications/TransferUniProtModifications/bin/Release/netcoreapp2.1/TransferUniProtModifications.dll",
        unixml=UNIPROTXML,
        protxml="data/combined.spritz.snpeff.protein.xml"
    output:
        protxml="data/combined.spritz.snpeff.protein.withmods.xml"
    shell:
        "dotnet {input.transfermods} -x {input.unixml} -y {input.protxml}"

rule transfer_modifications2:
    input:
        transfermods="TransferUniProtModifications/TransferUniProtModifications/bin/Release/netcoreapp2.1/TransferUniProtModifications.dll",
        unixml=UNIPROTXML,
        protxml="data/combined.spritz.isoform.protein.xml"
    output:
        protxml="data/combined.spritz.isoform.protein.withmods.xml"
    shell:
        "dotnet {input.transfermods} -x {input.unixml} -y {input.protxml}"

rule transfer_modifications3:
    input:
        transfermods="TransferUniProtModifications/TransferUniProtModifications/bin/Release/netcoreapp2.1/TransferUniProtModifications.dll",
        unixml=UNIPROTXML,
        protxml="data/GRCh38.86.protein.xml"
    output:
        protxml=
            "data/GRCh38.86.protein.withmods.xml",
    shell:
        "dotnet {input.transfermods} -x {input.unixml} -y {input.protxml}"

rule transfer_modifications4:
    input:
        transfermods="TransferUniProtModifications/TransferUniProtModifications/bin/Release/netcoreapp2.1/TransferUniProtModifications.dll",
        unixml=UNIPROTXML,
        protxml="data/combined.spritz.isoformvariants.protein.xml"
    output:
        protxml=
            "data/combined.spritz.isoformvariants.protein.withmods.xml"
    shell:
        "dotnet {input.transfermods} -x {input.unixml} -y {input.protxml}"

rule reference_protein_xml:
    """
    Create protein XML with sequences from the reference gene model.
    """
    input:
        "data/SnpEffDatabases.txt",
        snpeff="SnpEff/snpEff.jar",
        fa="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa",
    output:
        protxml="data/GRCh38.86.protein.xml"
    params:
        ref="GRCh38.86", # no isoform reconstruction
    resources:
        mem_mb=16000
    log:
        "data/GRCh38.86.spritz.log"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -nostats"
        " -xmlProt {output.protxml} {params.ref}) 2> {log}" # no isoforms, no variants

rule custom_protein_xml:
    """
    Create protein XML with sequences from the isoform discovery gene model.
    """
    input:
        "data/SnpEffDatabases.txt",
        snpeff="SnpEff/snpEff.jar",
        fa="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa",
        isoform_reconstruction="SnpEff/data/combined.sorted.filtered.withcds.gtf/genes.gtf"
    output:
        protxml="data/combined.spritz.isoform.protein.xml"
    params:
        ref="combined.sorted.filtered.withcds.gtf" # with isoforms
    resources:
        mem_mb=16000
    log:
        "data/combined.spritz.isoform.log"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -nostats"
        " -xmlProt {output.protxml} {params.ref}) 2> {log}" # isoforms, no variants
