rule unzip_ensembl:
    input:
        gfa="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
        gff="data/ensembl/Homo_sapiens.GRCh38.81.gff3.gz",
        pfa="data/ensembl/Homo_sapiens.GRCh38.pep.all.fa.gz",
        vcf="data/ensembl/common_all_20170710.vcf.gz",
        vcfidx="data/ensembl/common_all_20170710.vcf.idx.gz"
    output:
        gfa="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        gff="data/ensembl/Homo_sapiens.GRCh38.81.gff3",
        pfa="data/ensembl/Homo_sapiens.GRCh38.pep.all.fa",
        vcf="data/ensembl/common_all_20170710.vcf",
        vcfidx="data/ensembl/common_all_20170710.vcf.idx"
    shell:
        "gunzip -c {input.gfa} > {output.gfa} && "
        "gunzip -c {input.gff} > {output.gff} && "
        "gunzip -c {input.pfa} > {output.pfa} && "
        "gunzip -c {input.vcf} > {output.vcf} && "
        "gunzip -c {input.vcfidx} > {output.vcfidx} "

rule download_chromosome_mappings:
    output: "ChromosomeMappings/GRCh38_UCSC2ensembl.txt"
    shell: "git clone https://github.com/dpryan79/ChromosomeMappings.git"

rule reorder_genome_fasta:
    input: "data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output: "data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa"
    script: "../scripts/karyotypic_order.py"

rule convert_ucsc2ensembl:
    input:
        "data/ensembl/common_all_20170710.vcf",
        "ChromosomeMappings/GRCh38_UCSC2ensembl.txt"
    output:
        "data/ensembl/common_all_20170710.ensembl.vcf",
    script:
        "../scripts/convert_ucsc2ensembl.py"

rule index_ucsc2ensembl:
    input: "data/ensembl/common_all_20170710.ensembl.vcf"
    output: "data/ensembl/common_all_20170710.ensembl.vcf.idx"
    shell: "gatk IndexFeatureFile -F {input}"

rule filter_gff3:
    input: "data/ensembl/Homo_sapiens.GRCh38.81.gff3"
    output: "data/ensembl/202122.gff3"
    shell: "grep \"^#\|20\|^21\|^22\" \"data/ensembl/Homo_sapiens.GRCh38.81.gff3\" > \"data/ensembl/202122.gff3\""

rule filter_fa:
    input: "data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output: "data/ensembl/202122.fa"
    script: "../scripts/filter_fasta.py"

rule download_sras:
    output:
        temp("{dir}/{sra,[A-Z0-9]+}_1.fastq"), # constrain wildcards, so it doesn't soak up SRR######.trim_1.fastq
        temp("{dir}/{sra,[A-Z0-9]+}_2.fastq")
    log: "{dir}/{sra}.log"
    threads: 4
    shell:
        "fasterq-dump --progress --threads {threads} --split-files --outdir {wildcards.dir} {wildcards.sra} 2> {log}"

rule compress_fastqs:
    input:
        temp("{dir}/{sra,[A-Z0-9]+}_1.fastq"),
        temp("{dir}/{sra,[A-Z0-9]+}_2.fastq")
    output:
        "{dir}/{sra,[A-Z0-9]+}_1.fastq.gz",
        "{dir}/{sra,[A-Z0-9]+}_2.fastq.gz"
    shell:
        "gzip {input}"
