import pandas as pd
import glob
from itertools import chain
from itertools import product


def list_files(indir, extensions):
    extensions = extensions.split(",")
    files = [glob.glob(f"{indir}/*.{e}") for e in extensions]
    return list(chain.from_iterable(files))


def get_cov_seed(rep, sample2seed, general_seed):
    """
    Function that returns a seed for coverage
    By default each replicate has a unique seed
    If the param same_dist_rep is True, return a single seed for all replicates
    """
    try:
        same_dist_rep = config["same_dist_rep"]
        if same_dist_rep:
            return general_seed
    except KeyError:
        return sample2seed[int(rep)]


outdir = config["outdir"]

replicates = list(range(1, config["replicates"] + 1))

samples = []
if os.path.isfile(config["input_path"]):
    files = config["input_path"]
else:
    files = list_files(config["input_path"], "tsv")

for file in files:
    df = pd.read_csv(file, sep="\t")
    try:
        samples.append(list(set(df["sample"])))
    except KeyError:
        samples.append([os.path.basename(file).split(".")[0]])

samples = list(chain.from_iterable(samples))
samples_rep = list(product(samples, replicates))
renamed_samples = [
    "-".join([os.path.basename(sample[0]).split(".")[0], str(sample[1])])
    for sample in samples_rep
]


rule samples_table:
    input:
        f"{outdir}/cov.tsv",


def get_assembly_summary_path(wildcards):
    try:
        table = os.path.abspath(
            checkpoints.download_assemblies.get(**wildcards).output[0]
        )
    except AttributeError:
        table = os.path.abspath(config["asm_summary"])
    return table


rule get_samples_and_replicates:
    input:
        config["input_path"],
    output:
        temp(f"{outdir}/samples.tsv"),
    params:
        rep=config["replicates"],
        rep_sd=config["rep_sd"],
        seed=config["seed"],
    script:
        "get_rep_table.py"


checkpoint calculate_coverage:
    input:
        entry=f"{outdir}/samples.tsv",
        asm=lambda wildcards: get_assembly_summary_path(wildcards),
    output:
        temp(f"{outdir}/cov.tsv"),
    params:
        dist=config["dist"],
        mu=config["mu"],
        sigma=config["sigma"],
        total_bases=config["total_bases"],
        read_len=config["read_len"],
        pairing=config["paired"],
        rep_sd=config["rep_sd"],
        seed=config["seed"],
    log:
        f"logs/tables/cov.tsv",
    script:
        "calculate_cov.py"
