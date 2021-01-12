rule star:
    input: #TODO use expand of wildcards here reads_prefix?
        read1=reads[0],
        read2=reads[1]
    output:
        tx_align=temp("/".join([output_directory, samplename+"Aligned.toTranscriptome.out.bam"])),
        genome_align=temp("/".join([output_directory, samplename+"Aligned.out.bam"])),
        genome_align_sorted=temp("/".join([output_directory, samplename+"Aligned.sorted.bam"])),
        junc_file="/".join([output_directory, samplename+"SJ.out.tab"])
    params:
        index=config["resources"]["star_index"],
        gtf=config["resources"]["gtf"],
        prefix="/".join([output_directory, samplename]),
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
        --outSAMattrRGline ID:{samplename}  SM:{samplename}  PL:ILLUMINA \
        --outSAMattributes All \
        --outSAMunmapped Within \
        --outSAMprimaryFlag AllBestScore \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --chimOutType Junctions \
        --chimMainSegmentMultNmax 1 \
        --genomeLoad NoSharedMemory
        
        samtools sort -@ 15 -m 2G -o {output.genome_align_sorted} {output.genome_align} 
        # java -jar {params.picard} SortSam I={output.genome_align} \
        # O={output.genome_align_sorted} \
        # SO=coordinate
        """

rule rsem:
    input:
        rules.star.output.tx_align
    output:
        gene="/".join([output_directory, samplename+".genes.results"]),
        isoform="/".join([output_directory, samplename+".isoforms.results"])
    params:
        forward_prob="0",
        index=config["resources"]["rsem_index"],
        max_len="1000"
    threads: 15
    shell:
        """
        rsem-calculate-expression \
            --bam \
            --num-threads 15 \
            --fragment-length-max {params.max_len} \
            --no-bam-output \
            --paired-end \
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
        "/".join([output_directory, samplename+".idxstats"])
    threads: 1
    shell:
        """
        samtools idxstats {input} > {output}
        """