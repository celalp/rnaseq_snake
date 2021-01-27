rule intermediate_qc:
    input:
        rnaseq_metrics=rules.rnaseq_metrics.output.metrics,
        insert_size=rules.insert_size_metrics.output.metrics,
        alignment_metrics=rules.alignment_metrics.output,
        fastqc=rules.fastqc.output.fastqc_result,
        rnaseqc = rules.rnaseqc.output
    output:
        fake_qc=os.path.join(output_directory, "qc_metrics/{sampleID}qc_summary.txt")#a rule to put everything together
    group: "qc-metrics"
    threads: 2
    shell:
        """
        ls -ld {input.rnaseq_metrics} {input.insert_size} {input.alignment_metrics} \
        {input.fastqc} {input.rnaseqc} > {output.fake_qc}
        """

rule fastqc:
    input:
        reads=lambda wildcards: get_fastq(sample_df, wildcards.sampleID),
        qclibs=lambda wildcards: get_qc_libs(samples_df, wildcards.sampleID)
    threads: 4
    group: "qc-metrics"
    output:
        fastqc_result=expand(os.path.join(output_directory, "qc_metrics/{qclib}_fastqc/fastqc_data.txt"), qclib=qclibs)
    params:
        output_dir=os.path.join(output_directory, "qc_metrics")
    shell:
        "fastqc -t 4 -q --extract -o {params.output_dir} {input.reads} {input.read2}"

rule alignment_metrics:
    input:
        os.path.join(output_directory, "alignment/{sampleID}.mdup.bam")
    output: os.path.join(output_directory, "qc_metrics/{sampleID}.alignment_metrics.txt")
    params:
        picard=config["executables"]["picard"],
        fasta=config["resources"]["fasta"]
    threads: 4
    group: "qc-metrics"
    shell:
        """
        java -Xmx4096m -jar {params.picard} CollectAlignmentSummaryMetrics \
        I={input} R={params.fasta} O={output}        
        """

rule insert_size_metrics:
    input:
        os.path.join(output_directory, "alignment/{sampleID}.mdup.bam")
    output:
        metrics=os.path.join(output_directory, "qc_metrics/{sampleID}.insert_metrics.txt"),
        plot=os.path.join(output_directory, "qc_metrics/{sampleID}.insert_metrics.pdf")
    params:
        picard=config["executables"]["picard"]
    threads: 4
    group: "qc-metrics"
    shell:
        """
        java -Xmx4096m -jar {params.picard} CollectInsertSizeMetrics \
        I={input} O={output.metrics} H={output.plot}
        """

rule rnaseq_metrics:
    input:
        os.path.join(output_directory, "alignment/{sampleID}.mdup.bam")
    output:
        metrics=os.path.join(output_directory, "qc_metrics/{sampleID}.rnaseq_metrics.txt")
    params:
        ribo_int=config["resources"]["ribo_int"],
        refflat=config["resources"]["refflat"],
        strand="SECOND_READ_TRANSCRIPTION_STRAND",
        AS="true",
        picard=config["executables"]["picard"]
    threads: 4
    group: "qc-metrics"
    shell:
        """
        java -Xmx4096m -jar {params.picard} CollectRnaSeqMetrics \
        REF_FLAT={params.refflat} \
        RIBOSOMAL_INTERVALS={params.ribo_int} \
        STRAND_SPECIFICITY= {params.strand} \
        I={input} \
        O={output.metrics} AS={params.AS}
        """


rule rnaseqc:
    input: 
        bam = os.path.join(output_directory, "alignment/{sampleID}.mdup.bam"),
        collapsed_gtf = config["resources"]["collapsed_gtf"]
    output:
        os.path.join(output_directory, "qc_metrics/{sampleID}.metrics.tsv")
    params: 
        rnaseqc = config["executables"]["rnaseqc"],
        output_dir=os.path.join(output_directory, "qc_metrics")
    group: "qc-metrics"
    shell:
        """
        {params.rnaseqc} {input.collapsed_gtf} {input.bam} {params.output_dir} -s {wildcards.sampleID}
        rm -f {params.output_dir}/*gct
        """


