import os
import pandas as pd
import yaml

#load differential expression config yaml
DE = yaml.load(open('de_config.yaml'), Loader=yaml.BaseLoader)
    
de_output_files=[i.split(',')[1]+'vs'+i.split(',')[2] for i in DE['contrast']]
#['KO_35vsKO_14', 'KO_35vsWT_14', 'KO_35vsWT_35']

configfile: "config.yaml"

#create the main output directory
output_directory=config["root_dir"]

#load samples
sample_df = pd.read_csv(config["sample_file"], sep="\t")
sampleIDs=sample_df['sample'].tolist()

wildcard_constraints:
    sampleID = "|".join(str(v) for v in sample_df['sample'])

