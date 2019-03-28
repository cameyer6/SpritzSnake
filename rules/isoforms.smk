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
