rule download_protein_xml:
    output: "human.protein.xml.gz"
    shell:
        "python scripts/download_protein_xml.py | gzip -c > human.protein.xml.gz"

rule build_transfer_mods:
    output: "TransferUniProtModifications/TransferUniProtModifications/bin/Release/netcoreapp2.1/TransferUniProtModifications.dll"
    shell:
        """
        cd TransferUniProtModifications
        dotnet restore
        dotnet build -c Release TransferUniProtModifications.sln
        """

rule transfer_modifications:
    input:
        transfermods="TransferUniProtModifications/TransferUniProtModifications/bin/Release/netcoreapp2.1/TransferUniProtModifications.dll",
        unixml="human.protein.xml.gz",
        # unixml="TestData/Homo_sapiens_202022.xml.gz",
        protxml="TestData/{sample}_1.spritz.tr.snpeff.protein.xml",
    output:
        "TestData/{sample}_1.spritz.tr.snpeff.protein.withmods.xml",
        "TestData/{sample}_1.spritz.tr.snpeff.protein.fasta"
    shell:
        "dotnet {input.transfermods} -x {input.unixml} -y {input.protxml}"
