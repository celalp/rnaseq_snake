rule differential_expression:
    input:
        file_list=expand(os.path.join(output_directory, "quantification/{sampleID}"+quantification_output_extension), sampleID=sampleIDs),
        quantification_folder=os.path.join(output_directory, "quantification"),
        de_configfile=config['de_configfile'],
        samples_table=config['sample_file']
    params:
        output_dir="/".join([output_directory, 'diff_ex']),
    output:
        directory(os.path.join(output_directory, "diff_ex")
    shell:
        """
        Rscript scripts/differential_expression.R -s {input.samples_table} -y {input.de_configfile} \
        -l {input.quantification_folder} -d {params.output_dir}
        """
