import os
import pandas as pd

configfile: "config.yaml"

#samples=pd.read_csv(config["sample_file"], header=None, sep=" ")
#reads=samples.iloc[0,1].split(",")
#lib1=reads[0].split("/").pop().replace(".fastq.gz", "")
#lib2=reads[1].split("/").pop().replace(".fastq.gz", "")
#samplename = str(samples.iloc[0,0])
#tissue=samples.iloc[0,2]
#sample_description=samples.iloc[0,3].replace(",", "\t").replace(" ", "")

#create the main output directory
output_directory=config["root_dir"]

#load samples
sample_df = pd.read_csv(config["sample_file"], sep="\t")
sampleIDs=sample_df['sample'].tolist()

wildcard_constraints:
    sampleID = "|".join(str(v) for v in sample_df['sample'])

