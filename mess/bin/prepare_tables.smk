import random
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
random.seed(config["seed"])
replicates = list(range(1, config["replicates"] + 1))
samples = list_files(config["input_path"], "tsv")
samples_rep = list(product(samples, replicates))
renamed_samples = [
    "-".join([os.path.basename(sample[0]).split(".")[0], str(sample[1])])
    for sample in samples_rep
]
seed_list = random.sample(range(1, 1000000), len(renamed_samples))
seed_dict = dict(zip(renamed_samples, seed_list))


rule all_prepare_tables:
    input:
        expand("{outdir}/{sample}-cov.json", outdir=outdir, sample=renamed_samples),


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
        temp(f"{outdir}/{{sample}}.tsv"),
    params:
        outdir=outdir,
        rep=config["replicates"],
        sd_rep=config["sd_rep"],
        seed=config["seed"],
    script:
        """
        get_rep_tables.py
        """


rule calculate_coverage:
    input:
        entry=f"{outdir}/{{sample}}.tsv",
        asm=lambda wildcards: get_assembly_summary_path(wildcards),
    output:
        f"{outdir}/{{sample}}-cov.json",
    params:
        dist=config["dist"],
        mu=config["mu"],
        sigma=config["sigma"],
        total_bases=config["total_bases"],
        read_len=config["read_len"],
        pairing=config["paired"],
        rep_sd=config["rep_sd"],
        seed=lambda wildcards: get_cov_seed(wildcards.sample, seed_dict, config["seed"]),
    log:
        f"logs/tables/{{sample}}.tsv",
    shell:
        """
        calculate_cov.py -i {input.entry} -a {input.asm} \
        -d {params.dist} -m {params.mu} -s {params.sigma} \
        -n {params.total_bases} -l {params.read_len} -p {params.pairing} \
        -sd {params.rep_sd} --seed {params.seed} -o {output} 
        """
