rule differential_expression:
    input:
        #file_list=expand(os.path.join(output_directory, "quantification/{sampleID}.genes.results"), sampleID=sampleIDs),
        file_list=expand(os.path.join(output_directory, "quantification/{sampleID}"+DE['quantification_output_extension']), sampleID=sampleIDs),
        quantification_folder=os.path.join(output_directory, "quantification")      
    params:
        output_dir="/".join([output_directory, 'diff_ex']),
        de_extension=DE['quantification_output_extension'],
        sample_table=DE['sample_meta_table'],
        de_formula=DE['formula'],
        de_tooltype=DE['tool_type'],
        de_contrast='-'.join(DE['contrast'])
    output:
        expand(os.path.join(output_directory, "diff_ex/{comparison}.txt"), comparison=de_output_files)
    shell:
        """
        Rscript scripts/differential_expression.R -l {input.quantification_folder} -s {params.sample_table} \
        -t {params.de_tooltype} -f {params.de_formula} -d {params.output_dir} \
        -c {params.de_contrast} -e {params.de_extension}
        """