rule assemble_transcripts:
    input:
        bam="TestData/ERR315327_1.sorted.bam",
        gff="ensembl/202122.gff3"
    output: "TestData/ERR315327_1.sorted.gtf"
    threads: 12
    log: "TestData/ERR315327_1.sorted.gtf.log"
    shell: "stringtie {input.bam} -p {threads} -G {input.gff} -o {output} 2> {log}" # strandedness: --fr for forwared or --rf for reverse

# rule filter_gtf_entries_without_strand
#     input: "TestData/ERR315327_1.sorted.gtf"
#     output: "TestData/ERR315327_1.sorted.filtered.gtf"
#     script:

rule build_gtf_sharp:
    output: "GtfSharp/GtfSharp/bin/Release/netcoreapp2.1/GtfSharp.dll"
    shell:
        """
        cd GtfSharp
        dotnet restore
        dotnet build -c Release GtfSharp.sln
        """

rule filter_transcripts_add_cds:
    input:
        gtfsharp="GtfSharp/GtfSharp/bin/Release/netcoreapp2.1/GtfSharp.dll",
        gtf="TestData/ERR315327_1.sorted.gtf",
        fa="ensembl/202122.fa",
        refg="ensembl/202122.gff3"
    output:
        temp("TestData/ERR315327_1.sorted.filtered.gtf"),
        "TestData/ERR315327_1.sorted.filtered.withcds.gtf",
    shell:
        "dotnet {input.gtfsharp} -f {input.fa} -g {input.gtf} -r {input.refg}"
