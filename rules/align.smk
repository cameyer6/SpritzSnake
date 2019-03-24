rule directories:
    output: directory("ensembl/202122/")
    shell: "mkdir -p ensembl/202122/"

rule star_setup:
    output: "STAR-2.6.0c/bin/Linux_x86_64/STAR"
    shell:
        "wget https://github.com/alexdobin/STAR/archive/2.6.0c.tar.gz"
        "tar xvf 2.6.0c.tar.gz"
        "rm 2.6.0c.tar.gz"

rule star_genome_generate:
    input:
        star="STAR-2.7.0e/bin/Linux_x86_64/STAR",
        genomeDir=directory("ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic"),
        fa="ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa",
        gff="ensembl/202122.gff3"
    output:
        "ensembl/202122/SA"
    shell:
        "{input.star} --runMode genomeGenerate --runThreadN {threads} --genomeDir {input.genomeDir} "
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
