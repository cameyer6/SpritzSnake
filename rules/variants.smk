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
    input: "data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa"
    output: "data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa.fai"
    shell: "samtools faidx data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa"

rule dict_fa:
    input: "data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa"
    output: "data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.dict"
    shell: "gatk CreateSequenceDictionary -R {input} -O {output}"

rule tmpdir:
    output: temp(directory("tmp"))
    shell: "mkdir tmp"

rule hisat2_groupmark_bam:
    input:
        sorted="{dir}/combined.sorted.bam",
        tmp=directory("tmp")
    output:
        grouped=temp("{dir}/combined.sorted.grouped.bam"),
        groupedidx=temp("{dir}/combined.sorted.grouped.bam.bai"),
        marked="{dir}/combined.sorted.grouped.marked.bam",
        markedidx="{dir}/combined.sorted.grouped.marked.bam.bai",
        metrics="{dir}/combined.sorted.grouped.marked.metrics"
    resources:
        mem_mb=16000
    params:
        java_options="--java-options \"-Xmx{resources.mem_mb}M -Dsamjdk.compression_level=9\""
    shell:
        "gatk {params.java_options} AddOrReplaceReadGroups -PU platform  -PL illumina -SM sample -LB library -I {input.sorted} -O {output.grouped} -SO coordinate --TMP_DIR {input.tmp} && "
        "samtools index {output.grouped} && "
        "gatk {params.java_options} MarkDuplicates -I {output.grouped} -O {output.marked} -M {output.metrics} --TMP_DIR {input.tmp} -AS true &&"
        "samtools index {output.marked}"

# Checks if quality encoding is correct, and then splits n cigar reads
rule split_n_cigar_reads:
    input:
        bam="{dir}/combined.sorted.grouped.marked.bam",
        fa="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa",
        fai="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa.fai",
        fadict="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.dict",
        tmp=directory("tmp")
    output:
        fixed=temp("{dir}/combined.fixedQuals.bam"),
        split=temp("{dir}/combined.sorted.grouped.marked.split.bam"),
        splitidx=temp("{dir}/combined.sorted.grouped.marked.split.bam.bai")
    resources:
        mem_mb=16000
    params:
        java_options="--java-options \"-Xmx{resources.mem_mb}M -Dsamjdk.compression_level=9\""
    shell:
        "gatk {params.java_options} FixMisencodedBaseQualityReads -I {input.bam} -O {output.fixed} && "
        "gatk {params.java_options} SplitNCigarReads -R {input.fa} -I {output.fixed} -O {output.split} --tmp-dir {input.tmp} || " # fix and split
        "gatk {params.java_options} SplitNCigarReads -R {input.fa} -I {input.bam} -O {output.split} --tmp-dir {input.tmp}; " # or just split
        "samtools index {output.split}" # always index

rule base_recalibration:
    input:
        knownsites="data/ensembl/common_all_20170710.ensembl.vcf",
        knownsitesidx="data/ensembl/common_all_20170710.ensembl.vcf.idx",
        fa="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa",
        bam="{dir}/combined.sorted.grouped.marked.split.bam",
        tmp=directory("tmp")
    output:
        recaltable=temp("{dir}/combined.sorted.grouped.marked.split.recaltable"),
        recalbam=temp("{dir}/combined.sorted.grouped.marked.split.recal.bam")
    resources:
        mem_mb=16000
    params:
        java_options="--java-options \"-Xmx{resources.mem_mb}M -Dsamjdk.compression_level=9\""
    shell:
        """
        gatk {params.java_options} BaseRecalibrator -R {input.fa} -I {input.bam} --known-sites {input.knownsites} -O {output.recaltable} --tmp-dir {input.tmp}
        gatk {params.java_options} ApplyBQSR -R {input.fa} -I {input.bam} --bqsr-recal-file {output.recaltable} -O {output.recalbam} --tmp-dir {input.tmp}
        samtools index {output.recalbam}
        """

rule call_gvcf_varaints:
    input:
        knownsites="data/ensembl/common_all_20170710.ensembl.vcf",
        knownsitesidx="data/ensembl/common_all_20170710.ensembl.vcf.idx",
        fa="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa",
        bam="{dir}/combined.sorted.grouped.marked.split.recal.bam",
        tmp=directory("tmp")
    output: temp("{dir}/combined.sorted.grouped.marked.split.recal.g.vcf.gz"),
    threads: 12
    resources:
        mem_mb=16000
    params:
        java_options="--java-options \"-Xmx{resources.mem_mb}M -Dsamjdk.compression_level=9\""
    shell:
        "gatk {params.java_options} HaplotypeCaller"
        " --native-pair-hmm-threads {threads}"
        " -R {input.fa} -I {input.bam}"
        " --min-base-quality-score 20 --dont-use-soft-clipped-bases true"
        " --dbsnp {input.knownsites} -O {output} --tmp-dir {input.tmp}"
        " -ERC GVCF --max-mnp-distance 3 &&"
        "gatk IndexFeatureFile -F {output}"

rule call_vcf_variants:
    input:
        fa="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa",
        gvcf="{dir}/combined.sorted.grouped.marked.split.recal.g.vcf.gz",
        tmp=directory("tmp")
    output: "{dir}/combined.sorted.grouped.marked.split.recal.g.gt.vcf" # renamed in next rule
    resources:
        mem_mb=16000
    params:
        java_options="--java-options \"-Xmx{resources.mem_mb}M -Dsamjdk.compression_level=9\""
    shell:
        """
        gatk {params.java_options} GenotypeGVCFs -R {input.fa} -V {input.gvcf} -O {output} --tmp-dir {input.tmp}
        gatk IndexFeatureFile -F {output}
        """

rule final_vcf_naming:
    input: "{dir}/combined.sorted.grouped.marked.split.recal.g.gt.vcf"
    output: "{dir}/combined.spritz.vcf"
    shell: "mv {input} {output}"

rule filter_indels:
    input:
        fa="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa",
        vcf="{dir}/combined.spritz.vcf"
    output:
        "{dir}/combined.spritz.noindels.vcf"
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
    resources:
        mem_mb=16000
    shell:
        "java -Xmx{resources.mem_mb}M -jar {input.jar} databases > {output} && "
        "echo \"\n# {params.ref}\" >> SnpEff/snpEff.config && "
        "echo \"{params.ref}.genome : Human genome GRCh38 using RefSeq transcripts\" >> SnpEff/snpEff.config && "
        "echo \"{params.ref}.reference : ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/\" >> SnpEff/snpEff.config && "
        "echo \"\t{params.ref}.M.codonTable : Vertebrate_Mitochondrial\" >> SnpEff/snpEff.config && "
        "echo \"\t{params.ref}.MT.codonTable : Vertebrate_Mitochondrial\" >> SnpEff/snpEff.config"

rule variant_annotation_ref:
    input:
        "data/SnpEffDatabases.txt",
        snpeff="SnpEff/snpEff.jar",
        fa="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa",
        vcf="{dir}/combined.spritz.vcf",
    output:
        ann="{dir}/combined.spritz.snpeff.vcf",
        html="{dir}/combined.spritz.snpeff.html",
        genesummary="{dir}/combined.spritz.snpeff.genes.txt",
        protfa="{dir}/combined.spritz.snpeff.protein.fasta",
        protxml="{dir}/combined.spritz.snpeff.protein.xml"
    params:
        ref="GRCh38.86", # no isoform reconstruction
    resources:
        mem_mb=16000
    log:
        "{dir}/combined.spritz.snpeff.log"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -stats {output.html}"
        " -fastaProt {output.protfa} -xmlProt {output.protxml} "
        " {params.ref} {input.vcf}" # no isoforms, with variants
        " > {output.ann}) 2> {log}"

rule variant_annotation_custom:
    input:
        "data/SnpEffDatabases.txt",
        snpeff="SnpEff/snpEff.jar",
        fa="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa",
        vcf="data/combined.spritz.vcf",
        isoform_reconstruction="SnpEff/data/combined.sorted.filtered.withcds.gtf/genes.gtf"
    output:
        ann="data/combined.spritz.isoformvariants.vcf",
        html="data/combined.spritz.isoformvariants.html",
        genesummary="data/combined.spritz.isoformvariants.genes.txt",
        protfa="data/combined.spritz.isoformvariants.protein.fasta",
        protxml="data/combined.spritz.isoformvariants.protein.xml"
    params:
        ref="combined.sorted.filtered.withcds.gtf" # with isoforms
    resources:
        mem_mb=16000
    log:
        "data/combined.spritz.isoformvariants.log"
    shell:
        "(java -Xmx{resources.mem_mb}M -jar {input.snpeff} -v -stats {output.html}"
        " -fastaProt {output.protfa} -xmlProt {output.protxml}"
        " {params.ref} {input.vcf}" # with isoforms and variants
        " > {output.ann}) 2> {log}"

# rule cleanup_snpeff:
#     input:
#         "data/combined.spritz.snpeff.protein.xml",
#         "data/combined.spritz.isoform.protein.xml",
#         "data/GRCh38.86.protein.xml", # see proteogenomics.smk
#         "data/combined.spritz.isoformvariants.protein.xml" # see proteogenomics.smk
#     output:
#         temp("clean_snpeff")
#     shell:
#         # clean up SnpEff
#         "touch snpeff_complete && "
#         "cd SnpEff && git checkout -- snpEff.config && " # revert changes to config
#         "rm -r data/combined.sorted.filtered.withcds.gtf data/genomes/combined.sorted.filtered.withcds.gtf.fa && " # remove custom gene model files
#         "cd .. && rm data/SnpEffDatabases.txt" # remove database list to trigger the setup next time
