#!/usr/bin/env python

import sys, os, argparse
import subprocess
import yaml
from snakemake.utils import read_job_properties


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--depend", help="space seperated list of jobids to wait for")
parser.add_argument("positional", action="append", nargs="?")
args = parser.parse_args()


jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)

# default parameters defined in cluster_spec (accessed via snakemake read_job_properties)
cluster_param = job_properties["cluster"]
resources = job_properties["resources"]

if job_properties["type"] == "single":
    cluster_param["name"] = job_properties["rule"]
elif job_properties["type"] == "group":
    cluster_param["name"] = job_properties["groupid"]
    cluster_param["error"] = job_properties["groupid"]+'.err'
    cluster_param["output"] = job_properties["groupid"]+'.out'
else:
    raise NotImplementedError(
        "Don't know what to do with job_properties['type']=={}".format(job_properties['type'])
    )


if "threads" in job_properties:
    ppn = "ppn=" + str(job_properties["threads"])

# don't overwrite default parameters if defined in rule (or config file)
if ("threads" in job_properties) and ("threads" not in cluster_param):
    cluster_param["threads"] = job_properties["threads"]
if "mem" in resources: #and "mem" not in cluster_param:
     cluster_param["mem"] = str(resources["mem"]) + "g"

if "join" in cluster_param and cluster_param["join"]:
    cluster_param["join"] = "oe"
    # store job outputs from each rule to "pbs" sub-folder
    if "output" in cluster_param:
        if not os.path.isdir(cluster_param["output"]):
            os.mkdir(cluster_param["output"])#may need to edit the output file, now is snakemake folder

# qsub takes mem as mem,vmem
if "mem" in cluster_param:
    cluster_param["mem"] = "mem={},vmem={}".format(
        cluster_param["mem"], cluster_param["mem"]
    )

# check which system you are on and load command command_options
key_mapping_file = os.path.join(os.path.dirname(__file__), "key_mappings.yaml")
command_options = yaml.load(open(key_mapping_file), Loader=yaml.BaseLoader)
system = command_options["system"]
command = command_options[system]["command"]

key_mapping = command_options[system]["key_mapping"]

# construct command:
for key in key_mapping:
    if key in cluster_param:
        command += " "
        command += key_mapping[key].format(cluster_param[key])

depend = ""
if args.depend:
    jobids = ":".join(args.depend.split(" "))
    depend = ' -W "depend=afterok:{}" '.format(jobids)

if depend:
    command += depend

command += " {}".format(jobscript)

eprint("submit command: " + command)

try:
    res = subprocess.run(command, check=True, shell=True, stdout=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    raise e
res = res.stdout.decode()
jobid = res.strip()
print(jobid)