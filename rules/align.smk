import os

rule directories:
    output: directory("data/ensembl/202122/")
    shell: "mkdir -p data/ensembl/202122/"

rule star_setup:
    output: "STAR-2.6.0c/bin/Linux_x86_64/STAR"
    shell:
        "wget https://github.com/alexdobin/STAR/archive/2.6.0c.tar.gz"
        "tar xvf 2.6.0c.tar.gz"
        "rm 2.6.0c.tar.gz"

rule star_genome_generate:
    input:
        star="STAR-2.7.0e/bin/Linux_x86_64/STAR",
        genomeDir=directory("data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic"),
        fa="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa",
        gff="data/ensembl/202122.gff3"
    output:
        "data/ensembl/202122/SA"
    shell:
        "{input.star} --runMode genomeGenerate --runThreadN {threads} --genomeDir {input.genomeDir} "
        "--genomeFastaFiles {input.fa} --sjdbGTFfile {input.gff} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100"

rule hisat_genome:
    input:
        fa="data/ensembl/202122.fa",
        gtf="data/ensembl/202122.gff3"
    output: "data/ensembl/202122.1.ht2"
    shell: "hisat2-build data/ensembl/202122.fa data/ensembl/202122"

rule hisat2_splice_sites:
    input: "data/ensembl/202122.gff3"
    output: "data/ensembl/202122.splicesites.txt"
    shell: "hisat2_extract_splice_sites.py {input} > {output}"

def input_fq_args(fastqs):
    fqs=fastqs.split()
    if len(fqs) == 1:
        return f"-U {fqs[0]}"
    else:
        return f"-1 {fqs[0]} -2 {fqs[1]}"

rule hisat2_align_bam:
    input:
        "data/ensembl/202122.1.ht2",
        fq1="data/{sra}_1.fastq",
        fq2="data/{sra}_2.fastq",
        ss="data/ensembl/202122.splicesites.txt"
    output:
        sorted="data/{sra}.sorted.bam",
    threads: 12
    params:
        compression="9",
        tempprefix="data/{sra}.sorted"
    log: "data/{sra}.hisat2.log"
    shell:
        "(hisat2 -p {threads} -x data/ensembl/202122 -1 {input.fq1} -2 {input.fq2} --known-splicesite-infile {input.ss} | " # align the suckers
        "samtools view -h -F4 - | " # get mapped reads only
        "samtools sort -l {params.compression} -T {params.tempprefix} -o {output.sorted} -) 2> {log} && " # sort them
        "samtools index {output}"

rule hisat2_merge_bams:
    input:
        bams=expand("data/{sra}.sorted.bam", sra=config["sra"])
    output:
        sorted="data/combined.sorted.bam"
    params:
        compression="9",
        tempprefix="data/combined.sorted"
    log: "data/combined.sorted.log"
    shell:
        "(ls {input.bams} | "
        "{{ read firstbam; "
        "samtools view -h ""$firstbam""; "
        "while read bam; do samtools view ""$bam""; done; }} | "
        "samtools view -ubS - | "
        "samtools sort -l {params.compression} -T {params.tempprefix} -o {output.sorted} -) 2> {log} && "
        "samtools index {output.sorted}"
