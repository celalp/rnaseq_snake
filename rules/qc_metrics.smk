rule fastqc:
    input:
        read1=lambda wildcards: sample_df.loc[wildcards.sampleID, 'r1_path'],
        read2=lambda wildcards: sample_df.loc[wildcards.sampleID, 'r2_path']
    threads: 2
    group: "qc-metrics"
    output:
        report1=os.path.join(output_directory, "qc_metrics/{sampleID}_R1_fastqc/fastqc_data.txt"),
        report2=os.path.join(output_directory, "qc_metrics/{sampleID}_R2_fastqc/fastqc_data.txt")
    params:
        output_dir=os.path.join(output_directory, "qc_metrics")
    shell:
        "fastqc -t 1 -q --extract -o {params.output_dir} {input.read1} {input.read2}"

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
        rnaseqc = config["executables"]["rnaseqc"]
    group: "qc-metrics"
    shell:
        """
        {params.rnaseqc} {input.collapsed_gtf} {input.bam} {output_directory} -s {samplename}
        rm -f {output_directory}/*gct
        """
