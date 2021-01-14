rule star:
    input: #TODO use expand of wildcards here reads_prefix?
        read1=lambda wildcards: sample_df.loc[wildcards.sampleID, 'r1_path'],
        read2=lambda wildcards: sample_df.loc[wildcards.sampleID, 'r2_path']
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
    shell:
        """
        STAR --runMode alignReads \
        --runThreadN 15 \
        --readFilesCommand zcat \
        --readFilesIn {input.read1} {input.read2} \
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
        bam=temp(os.path.join(output_directory, "alignment/{sampleID}.mdup.bam")),
        metrics=temp(os.path.join(output_directory, "alignment/{sampleID}.duplicate.metrics"))
    params:
        picard=config["executables"]["picard"]
    threads: 10
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
    shell:
        """
        samtools idxstats {input} > {output.idx_stats}
        """