rule junctions:
    input:
        expr=ancient(rules.expression_outlier.output.gene),
        junctions=rules.star.output.junc_file,
        database=ancient(config["resources"]["gtex_db"]),
        txdb=config["resources"]["tx_db"],
        annotations=config["resources"]["junctions"],
        qc_check=ancient(rules.intermediate_qc.output.qc_sum)
        alia_file_rds=config["resources"]["alia_file"]
    params:
        tissue=tissue
    output:
        filtered = os.path.join(output_directory, samplename+".junctions.tsv"),
        alljuncs = os.path.join(output_directory, samplename + ".junctions.all.tsv")
    threads: 8
    shell:
        """
        Rscript scripts/junction_outlier.R -d {input.database} -t "{params.tissue}" -o {output.filtered} \
            -g {input.txdb} -a {input.annotations} -r {input.alia_file_rds} -j {input.junctions} -s {samplename} -e {input.expr} --ncores 8
        """