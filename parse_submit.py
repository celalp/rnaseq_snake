import pandas as pd
import os
import argparse as arg
import shutil
import yaml
import math

if __name__ == "__main__":
    parser = arg.ArgumentParser(description='Submit rna-seq alignment and quantitation jobs')
    parser.add_argument('-s', '--samples', type=str, help='tab separated samples file', action="store")
    parser.add_argument('-c', '--config', type=str, help='config.yaml file', action="store")
    parser.add_argument('-n', '--snakefile', type=str, help='snakefile', action="store")
    parser.add_argument('-b', '--bash_script', type=str, help='submission script', action="store")
    parser.add_argument('-o', '--outdir', type=str, help='name of the output directory', action="store")
    parser.add_argument('-r', '--resume', help='resubmit, resume after a failed job', action='store_true', default=False)
    args = parser.parse_args()

    samples = pd.read_csv(args.samples, header=None, sep="\t")
    if not os.path.isfile(args.config):
        raise FileNotFoundError("could not find the config file "+args.config)
    else:
        with open(args.config) as f:
            config = yaml.safe_load(f)
            f.close()

    # check sample files there's prob a better way of doing this
    for i in range(samples.shape[0]):
        read1 = samples.iloc[i, 1]
        read2 = samples.iloc[i, 2]
        if not os.path.isfile(read1):
            raise FileNotFoundError("could not find read1 for "+samples.iloc[i,0])
        if not os.path.isfile(str(read2)) and not pd.isnull(read2):
            raise FileNotFoundError("could not find read2 for "+samples.iloc[i,0])
        else:
            if pd.isnull(read2):
              print(samples.iloc[i,0] + " is a single end library")
            samplename=samples.iloc[i, 0]
            output_directory=args.outdir
            output_directory=os.path.abspath(output_directory)+"/"+samplename
            if os.path.isdir(output_directory) and not args.resume:
                print("Directory %s exists for same %s skipping" %(output_directory, samplename))
                continue
            else:
                os.makedirs(output_directory, exist_ok=True)
                config["samplename"] = samplename
                if pd.isnull(read2):
                    config["reads"] = {"read1": read1}
                else:
                    config["reads"] = {"read1": read1, "read2":read2}
                config["output_directory"] = output_directory

                # write the new config yaml file to output directory this is for reproducibility
                new_yaml = output_directory + "/config.yaml"
                with open(new_yaml, 'w') as outfile:
                    yaml.dump(config, outfile, default_flow_style=False)

                shutil.copy(args.snakefile, output_directory)
                command="qsub -v outdir={outdir} {script} -N {samplename} -e {outdir}/{samplename}.err -o {outdir}/{samplename}.out".\
                    format(samplename=samplename, script=args.bash_script, outdir=output_directory)

                print(command)
                os.system(command)




