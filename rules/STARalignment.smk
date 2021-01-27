rule star:
    input: 
        #read1=lambda wildcards: sample_df.loc[wildcards.sampleID, 'r1_path'],
        #read2=lambda wildcards: sample_df.loc[wildcards.sampleID, 'r2_path']
        reads=lambda wildcards: get_fastq(sample_df, wildcards.sampleID)
    output:
        tx_align=temp(os.path.join(output_directory, "alignment/{sampleID}Aligned.toTranscriptome.out.bam")),
        genome_align=temp(os.path.join(output_directory, "alignment/{sampleID}Aligned.out.bam")),
        genome_align_sorted=temp(os.path.join(output_directory, "alignment/{sampleID}Aligned.sorted.bam")),
        junc_file=os.path.join(output_directory, "alignment/{sampleID}SJ.out.tab")
    params:
        index=config["resources"]["star_index"],
        gtf=config["resources"]["gtf"],
        prefix=os.path.join(output_directory, "alignment/{sampleID}"),
        picard=config["executables"]["picard"]
    threads: 15
    group: "alignment-expression"
    shell:
        """
        STAR --runMode alignReads \
        --runThreadN 15 \
        --readFilesCommand zcat \
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
        --outSAMattrRGline ID:{wildcards.sampleID}  SM:{wildcards.sampleID}  PL:ILLUMINA \
        --outSAMattributes All \
        --outSAMunmapped Within \
        --outSAMprimaryFlag AllBestScore \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --chimOutType Junctions \
        --chimMainSegmentMultNmax 1 \
        --genomeLoad NoSharedMemory
        
        samtools sort -@ 15 -m 2G -o {output.genome_align_sorted} {output.genome_align} 
        """

rule mark_duplicates:
    input:
        rules.star.output.genome_align_sorted
    output:
        bam=os.path.join(output_directory, "alignment/{sampleID}.mdup.bam"),
        metrics=os.path.join(output_directory, "alignment/{sampleID}.duplicate.metrics")
    params:
        picard=config["executables"]["picard"]
    threads: 10
    group: "alignment-expression"
    shell:
        """
        java -Xmx24g -jar {params.picard} MarkDuplicates \
        I={input} O={output.bam} M={output.metrics} \
        ASSUME_SORT_ORDER=coordinate 
        
        samtools index -@ 10 {output.bam}
        """

rule samtools_idx:
    input:
        rules.mark_duplicates.output.bam
    output:
        idx_stats=os.path.join(output_directory, "alignment/{sampleID}.idxstats")
    threads: 1
    group: "alignment-expression"
    shell:
        """
        samtools idxstats {input} > {output.idx_stats}
        """

if config['option'].lower() == 'rsem':
    rule rsem:
        input:
            bam=rules.star.output.tx_align
            libtype=lambda wildcards: get_rsem_libtype(sample_df, wildcards.sampleID)
        output:
            gene=os.path.join(output_directory, "quantification/{sampleID}.genes.results"),
            isoform=os.path.join(output_directory, "quantification/{sampleID}.isoforms.results")#remove this? no use for D.E.
        params:
            forward_prob="0",
            index=config["resources"]["rsem_index"],
            max_len="1000",
            samples_dir=os.path.join(output_directory, "quantification/{sampleID}")
        threads: 15
        group: "alignment-expression"
        shell:
            """
            rsem-calculate-expression \
                --bam \
                --num-threads 15 \
                --fragment-length-max {params.max_len} \
                --no-bam-output \
                {input.libtype} \
                --estimate-rspd \
                --calc-ci \
                --forward-prob {params.forward_prob} \
                {input.bam} \
                {params.index} {params.samples_dir}
            """