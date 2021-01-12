rule fastqc:
    input:
        read1=reads[0],
        read2=reads[1]
    threads: 2
    output:
        report1="/".join([output_directory, lib1+"_fastqc/fastqc_data.txt"]),
        report2="/".join([output_directory, lib2+"_fastqc/fastqc_data.txt"])
    shell:
        "fastqc -t 1 -q --extract -o {output_directory} {input.read1} {input.read2}"

rule cal_q30:
    input: 
        read1 = reads[0],
        read2 = reads[1]
    output:
        qchist = os.path.join(output_directory, samplename + ".qchist.txt")
    threads: 1
    shell:
        """
        module load bbmap/37.33; reformat.sh -Xmx1g in1={input.read1} in2={input.read2} qchist={output.qchist}
        """

rule get_cov:
    input: 
        bam = rules.mark_duplicates.output.bam,
        genome_index = config["resources"]["fai"],
        bed = config["resources"]["exons_bed"]
    output:
        cov_file = os.path.join(output_directory, samplename + ".coverage.stat")
    threads: 8
    shell:
        """
        scripts/get_rnaseq_cov.sh {input.bam} {input.genome_index} {input.bed}
        """

rule alignment_metrics:
    input:
        rules.mark_duplicates.output.bam
    output: "/".join([output_directory, samplename+".alignment_metrics.txt"])
    params:
        picard=config["executables"]["picard"],
        fasta=config["resources"]["fasta"]
    threads: 4
    shell:
        """
        java -Xmx4096m -jar {params.picard} CollectAlignmentSummaryMetrics \
        I={input} R={params.fasta} O={output}        
        """

rule insert_size_metrics:
    input:
        rules.mark_duplicates.output.bam
    output:
        metrics="/".join([output_directory, samplename+".insert_metrics.txt"]),
        plot="/".join([output_directory, samplename+".insert_metrics.pdf"])
    params:
        picard=config["executables"]["picard"]
    threads: 4
    shell:
        """
        java -Xmx4096m -jar {params.picard} CollectInsertSizeMetrics \
        I={input} O={output.metrics} H={output.plot}
        """

rule rnaseq_metrics:
    input:
        rules.mark_duplicates.output.bam
    output:
        metrics="/".join([output_directory, samplename+".rnaseq_metrics.txt"]),
    params:
        ribo_int=config["resources"]["ribo_int"],
        refflat=config["resources"]["refflat"],
        strand="SECOND_READ_TRANSCRIPTION_STRAND",
        AS="true",
        picard=config["executables"]["picard"]
    threads: 4
    shell:
        """
        java -Xmx4096m -jar {params.picard} CollectRnaSeqMetrics \
        REF_FLAT={params.refflat} \
        RIBOSOMAL_INTERVALS={params.ribo_int} \
        STRAND_SPECIFICITY= {params.strand} \
        I={input} \
        O={output.metrics} AS={params.AS}
        """


rule ercc_sirv:
    input:
        gene=rules.rsem.output.gene,
        iso=rules.rsem.output.isoform,
        ercc_conc=config["resources"]["spikeins"]
    output:
        ercc="/".join([output_directory, samplename+".ercc.txt"]),
        sirv="/".join([output_directory, samplename+".sirv.txt"])
    threads: 2
    shell:
        """
        python3 scripts/spikeins.py --gene {input.gene} --iso {input.iso} --actual {input.ercc_conc} \
            --samplename {samplename} --analysis_path {output_directory}
        """

rule rnaseqc:
    input: 
        bam = rules.mark_duplicates.output.bam,
        collapsed_gtf = config["resources"]["collapsed_gtf"]
    output:
        os.path.join(output_directory, samplename +  ".metrics.tsv")
    params: 
        rnaseqc = config["executables"]["rnaseqc"]
    shell:
        """
        {params.rnaseqc} {input.collapsed_gtf} {input.bam} {output_directory} -s {samplename}
        rm -f {output_directory}/*gct
        """

rule intermediate_qc:
    input:
        ercc=rules.ercc_sirv.output.ercc,
        sirv=rules.ercc_sirv.output.sirv,
        rnaseq=rules.rnaseq_metrics.output.metrics,
        insert=rules.insert_size_metrics.output.metrics,
        alignment=rules.alignment_metrics.output,
        duplicates=rules.mark_duplicates.output.metrics,
        fastqc1=rules.fastqc.output.report1,
        fastqc2=rules.fastqc.output.report2,
        rmd=config["resources"]["rmd"],
        q30 = os.path.join(output_directory, samplename + ".qchist.txt"),
        coverage = ancient(os.path.join(output_directory, samplename + ".coverage.stat")),
        rnaseqc = rules.rnaseqc.output
    output:
        db="/".join([output_directory, samplename+".intermediate.db"]),
        qc_sum = os.path.join(output_directory, samplename + ".qcsum.txt")
#        report="/".join([output_directory, samplename+".intermediate_qc_report.html"])
    threads: 2
    shell:
        """
        python3 scripts/intermediate_qc.py --ercc {input.ercc} --sirv {input.sirv} \
            --rnaseq {input.rnaseq} --insert {input.insert} \
            --alignment {input.alignment} --duplicates {input.duplicates} \
            --fastqc1 {input.fastqc1} --fastqc2 {input.fastqc2} --analysis_path {output_directory} \
            --samplename {samplename}
        Rscript scripts/sum_qc_metrics.r -s {samplename} -d {output_directory}

#        Rscript scripts/generate_qc_report.R -s {samplename} -r {input.rmd} --analysis_path {output_directory}
        """