rule download_adapters:
    output:
        temp(directory("BBMap")),
        "adapters.fa"
    shell: "git clone --depth 1 https://github.com/BioInfoTools/BBMap.git && cp BBMap/resources/adapters.fa ."

rule skewer:
    input:
        fq1="{dir}/{sra}_1.fastq",
        fq2="{dir}/{sra}_2.fastq",
        adapters="adapters.fa"
    output:
        fq1="{dir}/trimmed/{sra}.trim_1.fastq",
        fq2="{dir}/trimmed/{sra}.trim_2.fastq"
    threads: 6
    log: "{dir}/trimmed/{sra}-trimmed.status"
    shell:
        "skewer -q 19 -o {dir}/trimmed/{wildcards.sra} -t {threads} -x {input.adapters} {input.fq1} {input.fq2} &> {log} && "
        "mv {dir}/trimmed/{wildcards.sra}-trimmed-pair1.fastq {dir}/trimmed/{wildcards.sra}.trim_1.fastq &&"
        "mv {dir}/trimmed/{wildcards.sra}-trimmed-pair2.fastq {dir}/trimmed/{wildcards.sra}.trim_2.fastq"

rule fastqc_analysis:
    input:
        fq1=["{dir}/{sra}_1.fastq", "{dir}/trimmed/{sra}.trim_1.fastq"],
        fq2=["{dir}/{sra}_2.fastq", "{dir}/trimmed/{sra}.trim_2.fastq"]
    output:
        fq1=["{dir}/{sra}_1_fastqc.html", "{dir}/{sra}_1_fastqc.zip", "{dir}/trimmed/{sra}.trim_1_fastqc.html", "{dir}/trimmed/{sra}.trim_1_fastqc.zip"],
        fq2=["{dir}/{sra}_2_fastqc.html", "{dir}/{sra}_2_fastqc.zip", "{dir}/trimmed/{sra}.trim_2_fastqc.html", "{dir}/trimmed/{sra}.trim_2_fastqc.zip"],
    log: "{dir}/{sra}.fastqc.log"
    threads: 6
    shell:
        "fastqc -t {threads} {input.fq1} {input.fq2} 2> {log}"
