rule all:
    input:
        "TestData/mapper0.bam"

rule download_genome_fasta:
    output:
        "ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    shell:
        "wget -O - ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | "
        "gunzip -c > ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

rule download_gene_model:
    output:
        "ensembl/Homo_sapiens.GRCh38.81.gff3"
    shell:
        "wget -O - ftp://ftp.ensembl.org/pub/release-81/gff3/homo_sapiens/Homo_sapiens.GRCh38.81.gff3.gz | "
        "gunzip -c > ensembl/Homo_sapiens.GRCh38.81.gff3"

rule download_protein_fasta:
    output:
        "ensembl/Homo_sapiens.GRCh38.pep.all.fa"
    shell:
        "wget -O - ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz | "
        "gunzip -c > ensembl/Homo_sapiens.GRCh38.pep.all.fa"

rule filter_gff3:
    input:
        "ensembl/Homo_sapiens.GRCh38.81.gff3"
    output:
        "ensembl/202122.gff3"
    shell:
        "grep \"^#\|20\|^21\|^22\" \"ensembl/Homo_sapiens.GRCh38.81.gff3\" > \"ensembl/202122.gff3\""

rule filter_fa:
    input:
        "ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output:
        "ensembl/202122.fa"
    script:
        "scripts/filter_fasta.py"


# rule download_sras:
#     input:
#         lambda wildcards: config["sra"][wildcards.sample]
#     output:
#         "fastqs/{sample}_1.fastq"
#         "fastqs/{sample}_2.fastq"
#     log:
#         "fastqs/{sample}.log"
#     threads: 4
#     shell:
#         "fasterq-dump --progress --threads {threads} --split-files --outdir fastqs 2> {log}"

rule directories:
    output: directory("ensembl/202122/")
    shell: "mkdir -p ensembl/202122/"

rule star_genome_generate:
    input:
        genomeDir="ensembl/202122/",
        fa="ensembl/202122.fa",
        gff="ensembl/202122.gff3"
    output:
        "ensembl/202122/SA"
    shell:
        "STAR-2.7.0e/bin/Linux_x86_64/STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {input.genomeDir} "
        "--genomeFastaFiles {input.fa} --sjdbGTFfile {input.gff} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100"

rule hisat_genome:
    input:
        fa="ensembl/202122.fa",
        gtf="ensembl/202122.gff3"
    output: "ensembl/202122.1.ht2"
    shell: "hisat2-build ensembl/202122.fa ensembl/202122"

# rule hisat2_align_bam:
#     input:
#         "ensembl/202122.1.ht2",
#         fq="TestData/mapper0.fastq"
#     output: "TestData/mapper0.bam"
#     shell: "hisat2 -x ensembl/202122 -q {input.fq} | samtools view -l 9 -b -o TestData/mapper0.bam"

rule hisat2_align_sam:
    input:
        "ensembl/202122.1.ht2",
        fq="TestData/mapper0.fastq",
        ss="ensembl/202122.splicesites.txt"
    output: "TestData/mapper0.sam"
    shell: "hisat2 -x ensembl/202122 -q {input.fq} --known-splicesite-infile {input.ss} > TestData/mapper0.sam"

rule hisat2_splice_sites:
    input: gff="ensembl/202122.gff3"
    output: ss="ensembl/202122.splicesites.txt"
    shell: "hisat2_extract_splice_sites.py {input.gff} > {output.ss}"
