rule download_adapters:
    output:
        temp(directory("BBMap")),
        "adapters.fa"
    shell: "git clone --depth 1 https://github.com/BioInfoTools/BBMap.git && cp BBMap/resources/adapters.fa ."

rule skewer:
    input:
        fq1="data/{sra}_1.fastq",
        fq2="data/{sra}_2.fastq",
        adapters="adapters.fa"
    output:
        fq1="data/trimmed/{sra}.trim_1.fastq",
        fq2="data/trimmed/{sra}.trim_2.fastq"
    threads: 6
    log: "data/trimmed/{sra}-trimmed.status"
    shell:
        "skewer -q 19 -o data/trimmed/{wildcards.sra} -t {threads} -x {input.adapters} {input.fq1} {input.fq2} &> {log} && "
        "mv data/trimmed/{wildcards.sra}-trimmed-pair1.fastq data/trimmed/{wildcards.sra}.trim_1.fastq &&"
        "mv data/trimmed/{wildcards.sra}-trimmed-pair2.fastq data/trimmed/{wildcards.sra}.trim_2.fastq"

rule fastqc_analysis:
    input:
        fq1=["data/{sra}_1.fastq", "data/trimmed/{sra}.trim_1.fastq"],
        fq2=["data/{sra}_2.fastq", "data/trimmed/{sra}.trim_2.fastq"]
    output:
        fq1=["data/{sra}_1_fastqc.html", "data/{sra}_1_fastqc.zip", "data/trimmed/{sra}.trim_1_fastqc.html", "data/trimmed/{sra}.trim_1_fastqc.zip"],
        fq2=["data/{sra}_2_fastqc.html", "data/{sra}_2_fastqc.zip", "data/trimmed/{sra}.trim_2_fastqc.html", "data/trimmed/{sra}.trim_2_fastqc.zip"],
    log: "data/{sra}.fastqc.log"
    threads: 6
    shell:
        "fastqc -t {threads} {input.fq1} {input.fq2} 2> {log}"
