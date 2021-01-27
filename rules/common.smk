import os
import pandas as pd

configfile: "config.yaml"

#create the main output directory
output_directory=config["root_dir"]

#different quantification tool output extensions
exten = {'rsem': '.genes.results', 'kallisto': 'abundance.h5', 
'salmon':'quant.sf.gz', 'sailfish': 'quant.sf'}
quantification_output_extension=exten[config['option'].lower()]

#load samples
sample_df = pd.read_csv(config["sample_file"], sep="\t")
sampleIDs=sample_df['sample'].tolist()

wildcard_constraints:
    sampleID = "|".join(str(v) for v in sample_df['sample'])

###helper functions

def get_fastq(samples_df, sampleid):
    """Get fastq files of given sample."""
    fastqs = samples_df.loc[sampleid, ['read1', 'read2']].dropna()
    if len(fastqs) == 2:
        return [fastqs.read1, fastqs.read2]
    elif len(fastqs) == 1:
        return [fastqs.read1]
    else:
        return "Missing fastq files"


def get_qc_libs(samples_df, sampleid):
    """Create qc_libs of given sample for qc analysis."""
    fastqs = samples_df.loc[sampleid, ['read1', 'read2']].dropna()
    lib1=fastqs.read1.split("/").pop().replace(".fastq.gz", "")
    if len(fastqs) == 2:
        lib2=fastqs.read2.split("/").pop().replace(".fastq.gz", "")
        return [lib1, lib2]
    elif len(fastqs) == 1:
        return [lib1]
    else:
        return "Missing fastq files"

def get_rsem_libtype(samples_df, sampleid):
    """For rsem input"""
    fastqs = samples_df.loc[sampleid, ['read1', 'read2']].dropna()
    if len(fastqs) == 2:
        libtype="--paired-end"
        return libtype
    elif len(fastqs) == 1:
        libtype=""
        return libtype
    else:
        return "Missing fastq files"