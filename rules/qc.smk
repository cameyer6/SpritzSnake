rule download_adapters:
    output:
        temp(directory("BBMap")),
        "data/qc/adapters.fa"
    shell: "git clone --depth 1 https://github.com/BioInfoTools/BBMap.git && cp BBMap/resources/adapters.fa data/qc"

rule skewer:
    input:
        fq1="{dir}/{sra}_1.fastq",
        fq2="{dir}/{sra}_2.fastq",
        adapters="data/qc/adapters.fa"
    output:
        fq1="{dir}/trimmed/{sra}.trim_1.fastq.gz",
        fq2="{dir}/trimmed/{sra}.trim_2.fastq.gz"
    threads: 12
    log: "{dir}/trimmed/{sra}-trimmed.status"
    params:
        quality=20
    shell:
        "skewer -q {params.quality} -o {wildcards.sra}/trimmed/{wildcards.sra}"
        " -t {threads} -x {input.adapters} {input.fq1} {input.fq2} &> {log} && "
        "mv {wildcards.dir}/trimmed/{wildcards.sra}-trimmed-pair1.fastq {wildcards.dir}/trimmed/{wildcards.sra}.trim_1.fastq &&"
        "mv {wildcards.dir}/trimmed/{wildcards.sra}-trimmed-pair2.fastq {wildcards.dir}/trimmed/{wildcards.sra}.trim_2.fastq &&"
        "gzip {wildcards.dir}/trimmed/{wildcards.sra}.trim_*.fastq"

rule fastqc_analysis:
    input:
        fq1=["{dir}/{sra}_1.fastq.gz", "{dir}/trimmed/{sra}.trim_1.fastq.gz"],
        fq2=["{dir}/{sra}_2.fastq.gz", "{dir}/trimmed/{sra}.trim_2.fastq.gz"]
    output:
        fq1=["{dir}/{sra}_1_fastqc.html", "{dir}/{sra}_1_fastqc.zip", "{dir}/trimmed/{sra}.trim_1_fastqc.html", "{dir}/trimmed/{sra}.trim_1_fastqc.zip"],
        fq2=["{dir}/{sra}_2_fastqc.html", "{dir}/{sra}_2_fastqc.zip", "{dir}/trimmed/{sra}.trim_2_fastqc.html", "{dir}/trimmed/{sra}.trim_2_fastqc.zip"],
    log: "{dir}/{sra}.fastqc.log"
    threads: 6
    shell:
        "fastqc -t {threads} {input.fq1} {input.fq2} 2> {log}"
