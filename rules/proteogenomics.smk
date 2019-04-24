rule download_protein_xml:
    output: "data/uniprot/human.protein.xml.gz"
    shell: "python scripts/download_protein_xml.py | gzip -c > {output}"

rule build_transfer_mods:
    output: "TransferUniProtModifications/TransferUniProtModifications/bin/Release/netcoreapp2.1/TransferUniProtModifications.dll"
    log: "data/TransferUniProtModifications.build.log"
    shell:
        "(cd TransferUniProtModifications && "
        "dotnet restore && "
        "dotnet build -c Release TransferUniProtModifications.sln) &> {log}"

rule transfer_modifications:
    input:
        transfermods="TransferUniProtModifications/TransferUniProtModifications/bin/Release/netcoreapp2.1/TransferUniProtModifications.dll",
        unixml="data/uniprot/human.protein.xml.gz",
        # unixml="data/Homo_sapiens_202022.xml.gz",
        protxml="data/combined.spritz.tr.snpeff.protein.xml",
    output:
        "data/combined.spritz.tr.snpeff.protein.withmods.xml",
        "data/combined.spritz.tr.snpeff.protein.fasta"
    shell:
        "dotnet {input.transfermods} -x {input.unixml} -y {input.protxml}"
