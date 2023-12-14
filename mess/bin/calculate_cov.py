#!/usr/bin/env python3

import pandas as pd
import numpy as np
import logging
import os
from humanfriendly import parse_size
import random


"""
Functions
"""


def get_lognormal_dist(df, mu, sigma):
    df["lognormal"] = np.random.lognormal(mean=mu, sigma=sigma, size=len(df))
    df["proportion"] = df["lognormal"] / df["lognormal"].sum()
    return df


def get_even_dist(df, cols):
    """
    function that calculates even abundances across taxonomy
    """
    abundances = df[cols].value_counts(normalize=True).reset_index()
    return df.merge(abundances)


"""
Main
"""
# Set logging
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%d %b %Y %H:%M:%S",
    filename=snakemake.log[0],
    level=logging.DEBUG,
)
# Set seed for distributions
np.random.seed(snakemake.params.seed)

# Set seed for read simulators
random.seed(snakemake.params.seed)

# Pairing
if snakemake.params.pairing:
    p = 2
else:
    p = 1


# Get table with assembly genomsizes and their taxonomy
entry_df = pd.read_csv(
    snakemake.input.entry, sep="\t", dtype={"entry": object, "chunk": object}
)
asm_df = pd.read_csv(snakemake.input.asm, sep="\t", dtype={"entry": object})
df = entry_df.merge(asm_df, on="entry")

# Calculate total base count
bases_per_sample = parse_size(snakemake.params.total_bases) // df["samplename"].count()

# Count number of entries
df["count"] = df.groupby(["samplename", "entry"])["entry"].transform("count")

# Calculate prportion with dist
if snakemake.params.dist == "even":
    df = get_even_dist(df, ["taxid"])
    df["bases"] = df["proportion"] * bases_per_sample
    df["reads"] = df["bases"] / (snakemake.params.read_len * p)
    df["cov_sim"] = df["bases"] / df["genome_size"]

elif snakemake.params.dist == "lognormal":
    df = get_lognormal_dist(df, snakemake.params.mu, snakemake.params.sigma)
    df["bases"] = df["proportion"] * bases_per_sample
    df["reads"] = df["bases"] / (snakemake.params.read_len * p)
    df["cov_sim"] = df["bases"] / df["genome_size"]

else:
    if "proportion" in entry_df.columns:
        df["proportion"] = df["proportion"] / df["count"]
        df["bases"] = df["proportion"] * bases_per_sample
        df["reads"] = df["bases"] / (snakemake.params.read_len * p)
        df["cov_sim"] = df["bases"] / df["genome_size"]

    elif "reads" in entry_df.columns:
        df["reads"] = df["reads"] / df["count"]
        df["bases"] = df["reads"] * (snakemake.params.read_len * p)
        bases_per_sample = df["bases"].sum()
        df["proportion"] = df["bases"] / bases_per_sample
        df["cov_sim"] = df["bases"] / df["genome_size"]

    elif "bases" in entry_df.columns:
        bases_per_sample = df["bases"].sum()
        df["bases"] = df["bases"] / df["count"]
        df["reads"] = df["bases"] / (snakemake.params.read_len * p)
        df["proportion"] = df["bases"] / bases_per_sample
        df["cov_sim"] = df["bases"] / df["genome_size"]

    elif "cov_sim" in entry_df.columns:
        df["cov_sim"] = df["cov_sim"] / df["count"]
        df["bases"] = df["cov_sim"] * df["genome_size"]
        df["reads"] = df["bases"] / (snakemake.params.read_len * p)
        df["proportion"] = df["bases"] / df["bases"].sum()


df["fasta"] = [os.path.basename(fa).split(".fna.gz")[0] for fa in df["path"]]
df["seed"] = random.sample(range(1, 1000000), len(df))

cols = [
    "samplename",
    "fasta",
    "genome_size",
    "contig_count",
    "bases",
    "reads",
    "proportion",
    "cov_sim",
    "chunk",
    "seed",
]
df = df.astype({"bases": int, "reads": int, "seed": int})
df[cols].to_csv(snakemake.log[0], sep="\t", index=None)  # type: ignore
df[cols].set_index(["samplename", "fasta"]).to_csv(snakemake.output[0], sep="\t")  # type: ignore
