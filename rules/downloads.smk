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

rule reorder_genome_fasta:
    output:
        "ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa"
    script:
        "../scripts/karyotypic_order.py"
