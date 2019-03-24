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

rule download_common_known_variants:
    output:
        "ensembl/common_all_20170710.vcf",
        "ensembl/common_all_20170710.vcf.idx"
    shell:
        "wget -O - ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/common_all_20170710.vcf.gz | "
        "gunzip -c > ensembl/common_all_20170710.vcf;"
        "gatk IndexFeatureFile -F ensembl/common_all_20170710.vcf"

rule download_chromosome_mappings:
    output:
        "ChromosomeMappings/GRCh38_ensembl2UCSC.txt"
    shell:
        "git clone https://github.com/dpryan79/ChromosomeMappings.git"

rule reorder_genome_fasta:
    input:
        "ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output:
        "ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa"
    script:
        "../scripts/karyotypic_order.py"

rule convert_ucsc2ensembl:
    input:
        "ensembl/common_all_20170710.vcf",
        "ChromosomeMappings/GRCh38_UCSC2ensembl.txt"
    output:
        "ensembl/common_all_20170710.ensembl.vcf",
    script:
        "../scripts/convert_ucsc2ensembl.py"

rule index_ucsc2ensembl:
    input: "ensembl/common_all_20170710.ensembl.vcf"
    output: "ensembl/common_all_20170710.ensembl.vcf.idx"
    shell: "gatk IndexFeatureFile -F {input}"

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
