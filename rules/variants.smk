rule download_snpeff:
    output: "SnpEff/snpEff.config", "SnpEff/snpEff.jar"
    log: "data/SnpEffInstall.log"
    shell:
        """
        (git clone --depth=1 https://github.com/smith-chem-wisc/SnpEff
        cd SnpEff
        mvn install:install-file -Dfile=lib/antlr-4.5.1-complete.jar -DgroupId=org.antlr -DartifactId=antlr -Dversion=4.5.1 -Dpackaging=jar
        mvn install:install-file -Dfile=lib/biojava3-core-3.0.7.jar -DgroupId=org.biojava -DartifactId=biojava3-core -Dversion=3.0.7 -Dpackaging=jar
        mvn install:install-file -Dfile=lib/biojava3-structure-3.0.7.jar -DgroupId=org.biojava -DartifactId=biojava3-structure -Dversion=3.0.7 -Dpackaging=jar
        export VERSION=4.3
        export VERSION_UND=`echo $VERSION | tr '.' '_'`
        mvn clean compile assembly:assembly
        mvn install:install-file -Dfile=target/SnpEff-$VERSION.jar -DgroupId=org.snpeff -DartifactId=SnpEff -Dversion=$VERSION -Dpackaging=jar -DgeneratePom=true --quiet
        cp target/SnpEff-$VERSION-jar-with-dependencies.jar snpEff.jar
        cd ..) &> {log}
        """

rule index_fa:
    input: "data/ensembl/202122.fa"
    output: "data/ensembl/202122.fa.fai"
    shell: "samtools faidx data/ensembl/202122.fa"

rule dict_fa:
    input: "data/ensembl/202122.fa"
    output: "data/ensembl/202122.dict"
    shell: "gatk CreateSequenceDictionary -R {input} -O {output}"

rule tmpdir:
    output: temp(directory("tmp"))
    shell: "mkdir tmp"

rule hisat2_group_bam:
    input:
        sorted="data/combined.sorted.bam",
        tmp=directory("tmp")
    output:
        grouped=temp("data/combined.sorted.grouped.bam"),
        groupedidx=temp("data/combined.sorted.grouped.bam.bai")
    shell:
        "gatk AddOrReplaceReadGroups -PU platform  -PL illumina -SM sample -LB library -I {input.sorted} -O {output.grouped} -SO coordinate --TMP_DIR tmp && "
        "samtools index {output.grouped}"

rule hisat2_mark_bam:
    input:
        grouped="data/combined.sorted.grouped.bam",
        tmp=directory("tmp")
    output:
        marked="data/combined.sorted.grouped.marked.bam",
        metrics="data/combined.sorted.grouped.marked.metrics"
    shell:
        "gatk MarkDuplicates -I {input.grouped} -O {output.marked} -M {output.metrics} --TMP_DIR tmp -AS true &&"
        "samtools index {output.marked}"

# Checks if quality encoding is correct, and then splits n cigar reads
rule split_n_cigar_reads:
    input:
        bam="data/combined.sorted.grouped.marked.bam",
        fa="data/ensembl/202122.fa",
        fai="data/ensembl/202122.fa.fai",
        fadict="data/ensembl/202122.dict"
    output:
        fixed=temp("data/combined.fixedQuals.bam"),
        split=temp("data/combined.sorted.grouped.marked.split.bam"),
        splitidx=temp("data/combined.sorted.grouped.marked.split.bam.bai")
    threads: 1
    shell:
        "gatk FixMisencodedBaseQualityReads -I {input.bam} -O {output.fixed} && "
        "gatk SplitNCigarReads -R {input.fa} -I {output.fixed} -O {output.split} || " # fix and split
        "gatk SplitNCigarReads -R {input.fa} -I {input.bam} -O {output.split}; " # or just split
        "samtools index {output.split}" # always index

rule base_recalibration:
    input:
        knownsites="data/ensembl/common_all_20170710.ensembl.vcf",
        knownsitesidx="data/ensembl/common_all_20170710.ensembl.vcf.idx",
        fa="data/ensembl/202122.fa",
        bam="data/combined.sorted.grouped.marked.split.bam"
    output:
        recaltable=temp("data/combined.sorted.grouped.marked.split.recaltable"),
        recalbam=temp("data/combined.sorted.grouped.marked.split.recal.bam")
    threads: 1
    shell:
        """
        gatk BaseRecalibrator -R {input.fa} -I {input.bam} --known-sites {input.knownsites} -O {output.recaltable}
        gatk ApplyBQSR -R {input.fa} -I {input.bam} --bqsr-recal-file {output.recaltable} -O {output.recalbam}
        samtools index {output.recalbam}
        """

rule call_gvcf_varaints:
    input:
        knownsites="data/ensembl/common_all_20170710.ensembl.vcf",
        knownsitesidx="data/ensembl/common_all_20170710.ensembl.vcf.idx",
        fa="data/ensembl/202122.fa",
        bam="data/combined.sorted.grouped.marked.split.recal.bam"
    output: temp("data/combined.sorted.grouped.marked.split.recal.g.vcf.gz"),
    threads: 4
    shell:
        "gatk HaplotypeCaller"
        " --native-pair-hmm-threads {threads}"
        " -R {input.fa} -I {input.bam}"
        " --min-base-quality-score 20 --dont-use-soft-clipped-bases true"
        " --dbsnp {input.knownsites} -O {output}"
        " -ERC GVCF --max-mnp-distance 3 &&"
        "gatk IndexFeatureFile -F {output}"

rule call_vcf_variants:
    input:
        fa="data/ensembl/202122.fa",
        gvcf="data/combined.sorted.grouped.marked.split.recal.g.vcf.gz",
    output: "data/combined.sorted.grouped.marked.split.recal.g.gt.vcf" # renamed in next rule
    shell:
        """
        gatk GenotypeGVCFs -R {input.fa} -V {input.gvcf} -O {output}
        gatk IndexFeatureFile -F {output}
        """

rule final_vcf_naming:
    input: "data/combined.sorted.grouped.marked.split.recal.g.gt.vcf"
    output: "data/combined.spritz.vcf"
    shell: "mv {input} {output}"

rule filter_indels:
    input:
        fa="data/ensembl/202122.fa",
        vcf="data/combined.spritz.vcf"
    output:
        "data/combined.spritz.noindels.vcf"
    shell:
        """
        gatk SelectVariants --select-type-to-exclude INDEL -R {input.fa} -V {input.vcf} -O {output}
        gatk IndexFeatureFile -F {output}
        """

# rule snpeff_data_folder:
#     output: directory("SnpEff/data")
#     shell: "mkdir SnpEff/data"

rule snpeff_database_setup:
    input:
        # dir="SnpEff/data",
        jar="SnpEff/snpEff.jar",
        config="SnpEff/snpEff.config"
    output:
        "data/SnpEffDatabases.txt"
    params:
        ref="GRCh38.86"
    shell:
        "java -Xmx2000M -jar {input.jar} databases > {output} &&"
        "echo \"\n# {params.ref}\" >> SnpEff/snpEff.config"
        "echo \"# {params.ref}.genome : Human genome GRCh38 using RefSeq transcripts\" >> SnpEff/snpEff.config"
        "echo \"# {params.ref}.reference : ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/\" >> SnpEff/snpEff.config"
        "echo \"# {params.ref}.M.codonTable : Vertebrate_Mitochondrial\" >> SnpEff/snpEff.config"
        "echo \"# {params.ref}.MT.codonTable : Vertebrate_Mitochondrial\" >> SnpEff/snpEff.config"

rule variant_annotation_ref:
    input:
        "data/SnpEffDatabases.txt",
        snpeff="SnpEff/snpEff.jar",
        fa="data/ensembl/202122.fa",
        vcf="data/combined.spritz.vcf",
    output:
        ann="data/combined.spritz.snpeff.vcf",
        html="data/combined.spritz.snpeff.html",
        genesummary="data/combined.spritz.snpeff.genes.txt",
        protfa="data/combined.spritz.snpeff.protein.fasta",
        protxml="data/combined.spritz.snpeff.protein.xml",
    params:
        ref="GRCh38.86"
    log:
        "data/combined.spritz.snpeff.log"
    shell:
        "(java -Xmx5000M -jar {input.snpeff} -v -stats {output.html}"
        " -fastaProt {output.protfa} -xmlProt {output.protxml} {params.ref}"
        " {input.vcf} > {output.ann}) 2> {log}"

rule variant_annotation_custom:
    input:
        "data/SnpEffDatabases.txt",
        snpeff="SnpEff/snpEff.jar",
        fa="data/ensembl/202122.fa",
        vcf="data/combined.spritz.vcf",
        trigger_isoform_reconstruction="SnpEff/data/combined.sorted.filtered.withcds.gtf/genes.gtf"
    output:
        ann="data/combined.spritz.tr.snpeff.vcf",
        html="data/combined.spritz.tr.snpeff.html",
        genesummary="data/combined.spritz.tr.snpeff.genes.txt",
        protfa="data/combined.spritz.tr.snpeff.protein.fasta",
        protxml="data/combined.spritz.tr.snpeff.protein.xml",
    params:
        ref="GRCh38.86"
    log:
        "data/combined.spritz.tr.snpeff.log"
    shell:
        "(java -Xmx5000M -jar {input.snpeff} -v -stats {output.html}"
        " -fastaProt {output.protfa} -xmlProt {output.protxml} {params.ref}"
        " {input.vcf} > {output.ann}) 2> {log}"
