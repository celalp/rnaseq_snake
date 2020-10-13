import os

configfile: "config.yaml"

reads1=config["reads"]["read1"]
lib1=reads1.split("/").pop().replace(".fastq.gz", "")

if "read2" in config["reads"].keys():
    reads2=config["reads"]["read2"]
    lib2=reads2.split("/").pop().replace(".fastq.gz", "")
    libtype="--paired-end"
    reads=[reads1, reads2]
    libs=[lib1, lib2]
else:
    libtype=""
    reads=reads1
    libs=[lib1]

if reads1.endswith(".gz"):
    compressed="--readFilesCommand zcat"
else:
    compressed=""

samplename=config["samplename"]

output_directory=config["output_directory"]


rule done:
    input:
        fastqc=expand("{output_directory}/{libname}_fastqc/fastqc_data.txt", output_directory=output_directory, libname=libs),
        alignment="/".join([output_directory, samplename+".alignment_metrics.txt"]),
        inserts="/".join([output_directory, samplename+".insert_metrics.txt"]),
        duplicates="/".join([output_directory, samplename+".duplicate.metrics"]),
        rnaseq="/".join([output_directory, samplename+".rnaseq_metrics.txt"]),
        bam="/".join([output_directory, samplename+".mdup.bam"]),
        gene="/".join([output_directory, samplename + ".genes.results"]),
        isoform = "/".join([output_directory, samplename + ".isoforms.results"])
    shell:
       """
       rm -rf {samplename}Log.out {samplename}Log.final.out {samplename}Log.progress.out \
          {samplename}_STARgenome {samplename}_STARpass1 {samplename}_STARtmp
       """

rule star:
    input:
        reads=reads
    output:
        tx_align = "/".join([output_directory, samplename + "Aligned.toTranscriptome.out.bam"]),
        genome_align = temp("/".join([output_directory, samplename + "Aligned.out.bam"])),
        genome_align_sorted = "/".join([output_directory, samplename + "Aligned.sorted.bam"]),
        junc_file = "/".join([output_directory, samplename + "SJ.out.tab"])
    params:
        index = config["resources"]["star_index"],
        gtf = config["resources"]["gtf"],
        prefix = "/".join([output_directory, samplename]),
        picard = config["executables"]["picard"]
    threads: 10
    shell:
        """
        STAR --runMode alignReads \
        --runThreadN 10 \
        {compressed} \
        --readFilesIn {input.reads} \
        --genomeDir {params.index} \
        --outFileNamePrefix {params.prefix} \
        --twopassMode Basic \
        --sjdbGTFfile {params.gtf} \
        --outFilterType BySJout \
        --limitSjdbInsertNsj 1200000 \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs None \
        --alignSoftClipAtReferenceEnds Yes \
        --quantMode TranscriptomeSAM GeneCounts \
        --outSAMtype BAM Unsorted \
        --outSAMattrRGline ID:{samplename}  SM:{samplename}  PL:ILLUMINA \
        --outSAMattributes All \
        --outSAMunmapped Within \
        --outSAMprimaryFlag AllBestScore \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --chimOutType Junctions \
        --chimMainSegmentMultNmax 1 \
        --genomeLoad NoSharedMemory

        java -jar {params.picard} SortSam I={output.genome_align} \
        O={output.genome_align_sorted} \
        SO=coordinate
        """

rule fastqc:
    input:
        reads=reads
    output:
        report1="/".join([output_directory, lib1+"_fastqc/fastqc_data.txt"]),
        report2="/".join([output_directory, lib2+"_fastqc/fastqc_data.txt"])
    threads: 10
    shell:
        "fastqc -t 10 -q --extract -o {output_directory} {input.reads}"

rule rsem:
    input:
        rules.star.output.tx_align
    output:
        gene = "/".join([output_directory, samplename + ".genes.results"]),
        isoform = "/".join([output_directory, samplename + ".isoforms.results"])
    params:
        forward_prob = "0",
        index = config["resources"]["rsem_index"],
        max_len = "1000",
    threads: 10
    shell:
        """
        rsem-calculate-expression \
            --bam \
            --num-threads 10 \
            --fragment-length-max {params.max_len} \
            --no-bam-output \
            {libtype} \
            --estimate-rspd \
            --calc-ci \
            --forward-prob {params.forward_prob} \
            {input} \
            {params.index} {output_directory}/{samplename}
        """


rule mark_duplicates:
    input:
        rules.star.output.genome_align_sorted
    output:
        bam="/".join([output_directory, samplename+".mdup.bam"]),
        metrics="/".join([output_directory, samplename+".duplicate.metrics"])
    params:
        picard=config["executables"]["picard"]
    threads: 2
    shell:
        """
        java -Xmx4096m -jar {params.picard} MarkDuplicates \
        I={input} O={output.bam} M={output.metrics} \
        ASSUME_SORT_ORDER=coordinate
        
        samtools index -@ 10 {output.bam}
        """


rule alignment_metrics:
    input:
        rules.mark_duplicates.output.bam
    output: "/".join([output_directory, samplename+".alignment_metrics.txt"])
    params:
        picard=config["executables"]["picard"],
        fasta=config["resources"]["fasta"]
    threads: 2
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
    threads: 2
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
        refflat=config["resources"]["refflat"],
        strand="SECOND_READ_TRANSCRIPTION_STRAND",
        AS="true",
        picard=config["executables"]["picard"]
    threads: 2
    shell:
        """
        java -Xmx4096m -jar {params.picard} CollectRnaSeqMetrics \
        REF_FLAT={params.refflat} \
        STRAND_SPECIFICITY= {params.strand} \
        I={input} \
        O={output.metrics} AS={params.AS}
        """

