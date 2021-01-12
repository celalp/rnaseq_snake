rule expression_outlier:
    input:
        gene_expressions=rules.rsem.output.gene,
        isoform_expressions=rules.rsem.output.isoform,
        database=ancient(config["resources"]["gtex_db"])
    params:
        tissue=tissue,
        output="/".join([output_directory, samplename])
    output:
        gene="/".join([output_directory, samplename+".gene.tsv"]),
        tx="/".join([output_directory, samplename+".tx.tsv"])
    threads: 5
    shell:
        """
        python3 scripts/expression_outlier.py -g {input.gene_expressions} \
            -i {input.isoform_expressions} -d {input.database} \
            -t "{params.tissue}" -o {params.output}
        """