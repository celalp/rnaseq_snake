rule split_cigar_reads:
    input:
        bam=rules.mark_duplicates.output.bam,
        qc_check=ancient(rules.intermediate_qc.output.qc_sum),
        junctions=rules.junctions.output.filtered
    output:
        temp("/".join([output_directory, samplename+".split.bam"]))
    params:
        fasta=config["resources"]["fasta"],
        gatk=config["executables"]["gatk"]
    threads: 10
    shell:
        """
        java -Xmx40960m -jar {params.gatk} \
        -T SplitNCigarReads \
        -R {params.fasta} \
        -I {input.bam} \
        -o {output} \
        -rf ReassignOneMappingQuality \
        -RMQF 255 \
        -RMQT 60 \
        --read_filter BadCigar --read_filter NotPrimaryAlignment \
        --logging_level ERROR \
        -U ALLOW_N_CIGAR_READS
        
        samtools index {output}
        """

rule recalibrate:
    input: rules.split_cigar_reads.output
    output:
        step1="/".join([output_directory, samplename+".recal_data.table"]),
        step2="/".join([output_directory, samplename+".post_recal_data.table"])
    params:
        gatk=config["executables"]["gatk"],
        fasta=config["resources"]["fasta"],
        sites1=config["resources"]["dbsnp"],
        sites2=config["resources"]["indels"]
    threads: 16
    shell:
        """
        java -Xmx40960m -jar {params.gatk} \
        -T BaseRecalibrator \
        -R {params.fasta} \
        -I {input} \
        -knownSites {params.sites1} \
        -knownSites {params.sites2} \
        --logging_level ERROR \
        -o {output.step1} \
        -nct 16
    
        java -Xmx40960m -jar {params.gatk} \
        -T BaseRecalibrator \
        -R {params.fasta} \
        -I {input} \
        -knownSites {params.sites1} \
        -knownSites {params.sites2} \
        -BQSR {output.step1} \
        --logging_level ERROR \
        -o {output.step2} \
        -nct 16
        """

rule print_reads:
    input:
        bam=rules.split_cigar_reads.output,
        recalibrate=rules.recalibrate.output.step1
    output:
        temp("/".join([output_directory, samplename+".recal_reads.bam"]))
    params:
        gatk=config["executables"]["gatk"],
        fasta=config["resources"]["fasta"]
    threads: 16
    shell:
        """
        java -Xmx64g -jar {params.gatk} \
        -T PrintReads \
        -R {params.fasta} \
        -I {input.bam} \
        -BQSR {input.recalibrate} \
        --logging_level ERROR \
        -o {output} \
        -nct 8
        """


rule haplotype_caller:
    input: rules.print_reads.output
    output: "/".join([output_directory, samplename+".vcf.gz"])
    params:
       gatk=config["executables"]["gatk"],
       fasta=config["resources"]["fasta"]
    threads: 16
    shell:
        """
         java -Xmx40960m -jar {params.gatk} \
        -T HaplotypeCaller \
        -R {params.fasta} \
        -I {input} \
        --annotation MappingQualityRankSumTest \
        --annotation MappingQualityZero \
        --annotation QualByDepth \
        --annotation ReadPosRankSumTest \
        --annotation RMSMappingQuality \
        --annotation BaseQualityRankSumTest \
        --annotation FisherStrand \
        --annotation GCContent \
        --annotation HaplotypeScore \
        --annotation HomopolymerRun \
        --annotation DepthPerAlleleBySample \
        --annotation Coverage \
        --annotation ClippingRankSumTest \
        --annotation DepthPerSampleHC \
        -dontUseSoftClippedBases \
        -stand_call_conf 20.0 \
        --logging_level ERROR \
        -o {output} \
        -nct 16
        """

rule variant_filter:
    input:
        vcf=rules.haplotype_caller.output,
        #filter_bed=rules.prepare_filter_bed.output,
        rna_edit_bed=config["resources"]["rna_edit"]
    output:
        qual_filt="/".join([output_directory, samplename+".qual_filt.vcf.gz"]),
        no_rna_edit=temp("/".join([output_directory, samplename+".no_rnaedit.vcf.gz"])),
        decompose="/".join([output_directory, samplename+".filtered.vcf"])
    params:
        gatk=config["executables"]["gatk"],
        fasta=config["resources"]["fasta"]
    threads: 16
    shell:
        """
        java -Xmx40960m -jar {params.gatk} \
        -T VariantFiltration \
        -R {params.fasta} \
        -V {input.vcf} \
        -window 35 \
        -cluster 3 \
        -filterName FS \
        -filter "FS > 30.0" \
        -filterName QD \
        -filter "QD < 2.0" \
        -o {output.qual_filt} 
        
        bcftools view -T {input.rna_edit_bed} -f "PASS" -o {output.no_rna_edit} -O z {output.qual_filt}
        
        gunzip -c {output.no_rna_edit} | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -s - | awk '{{ gsub("./-65", "./."); print $0 }}'| bgzip -c > {output.decompose}
        """

rule vep:
    input: rules.variant_filter.output.decompose
    output: "/".join([output_directory, samplename+".annotated.vcf.gz"])
    params:
        vep_cache=config["resources"]["vep_cache"],
        vep_ver=config["resources"]["vep_ver"],
        fasta=config["resources"]["fasta"],
        plugin_dir=config["resources"]["vep_plugin_dir"]
    threads: 16
    shell:
        """
        vep -o {output} -i {input} --species homo_sapiens   \
        --cache --offline --dir {params.vep_cache} --cache_version {params.vep_ver} \
        --dir_plugins {params.plugin_dir} --fork 16 --stats_text --af_exac \
        --sift b --polyphen b --shift_hgvs 1 --terms SO --vcf --everything \
        --fasta {params.fasta} --clin_sig_allele  --buffer_size 500 --compress_output bgzip
        
        tabix {output}
        """

rule pedfile:
    output:
        temp("/".join([output_directory, samplename+".ped"]))
    threads: 1
    shell:
        """
        echo -e "{sample_description}" > {output}
        """

rule vcf2db:
    input:
        vcf=rules.vep.output,
        ped=rules.pedfile.output
    output:
        "/".join([output_directory, samplename+".gemini.db"])
    threads: 5
    shell:
          """
          vcf2db.py {input.vcf} {input.ped} {output} --expand gt_types --expand gt_ref_depths \
          --expand gt_alt_depths --expand gt_quals \
          --expand gt_depths --expand gt_alt_freqs 
          """

rule bcftools:
    input:
        rules.variant_filter.output.decompose
    output:
        "/".join([output_directory, samplename+".bcftools.stats"])
    threads: 1
    shell:
        """
        bcftools stats {input} > {output}
        """

