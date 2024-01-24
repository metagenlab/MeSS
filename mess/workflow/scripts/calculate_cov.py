import pandas as pd
import numpy as np
from humanfriendly import parse_size
import random


"""
Functions
"""


def get_lognormal_dist(df, mu, sigma):
    """
    function for simulating lognormal sequence distribution
    """
    df["lognormal"] = np.random.lognormal(mean=mu, sigma=sigma, size=len(df))
    df["proportion"] = df["lognormal"] / df["lognormal"].sum()
    return df


def get_even_dist(df, cols):
    """
    function that calculates even abundances across taxonomy
    """
    abundances = (
        df[cols]
        .value_counts(normalize=True)
        .reset_index()
        .rename(columns={"proportion": "abundance"}, inplace=True)
    )
    return df.merge(abundances)


"""
Main
"""
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
entry_df = pd.read_csv(snakemake.input.df, sep="\t")

asm_df = pd.read_csv(snakemake.input.asm, sep="\t")
same_cols = list(np.intersect1d(entry_df.columns, asm_df.columns))
df = pd.merge(entry_df, asm_df, how="left", on=same_cols)


# Get base count per sample
bases = parse_size(snakemake.params.bases)

# Count number of entries
df["count"] = df.groupby(["samplename", "entry"])["entry"].transform("count")

# Calculate prportion with dist
if snakemake.params.dist == "even":
    df = get_even_dist(df, ["taxid"])
    df["sum_genome_size"] = df.groupby("samplename")["genome_size"].transform("sum")
    df["sum_cov"] = bases / df["sum_genome_size"]
    df["cov_sim"] = df["sum_cov"] * df["abundance"]


elif snakemake.params.dist == "lognormal":
    df = get_lognormal_dist(df, snakemake.params.mu, snakemake.params.sigma)
    df["bases"] = df["proportion"] * bases
    df["reads"] = df["bases"] / (snakemake.params.read_len * p)
    df["cov_sim"] = df["bases"] / df["genome_size"]
    df["sum_cov"] = df.groupby("samplename")["cov_sim"].transform("sum")
    df["abundance"] = df["cov_sim"] / df["sum_cov"]
    df["sum_bases"] = df.groupby("samplename")["bases"].transform("sum")
    df["reads"] = df["bases"] / (snakemake.params.read_len * p)
    df["proportion"] = df["bases"] / df["sum_bases"]
else:
    if "abundance" in entry_df.columns:
        df["sum_genome_size"] = df.groupby("samplename")["genome_size"].transform("sum")
        df["sum_cov"] = bases / df["sum_genome_size"]
        df["cov_sim"] = df["sum_cov"] * df["abundance"]
        df["bases"] = df["cov_sim"] * df["genome_size"]
        df["proportion"] = df["bases"] / df["sum_bases"]
        df["reads"] = df["bases"] / (snakemake.params.read_len * p)

    if "proportion" in entry_df.columns:
        df["proportion"] = df["proportion"] / df["count"]
        df["bases"] = df["proportion"] * bases
        df["reads"] = df["bases"] / (snakemake.params.read_len * p)
        df["cov_sim"] = df["bases"] / df["genome_size"]
        df["sum_cov"] = df.groupby("samplename")["cov_sim"].transform("sum")
        df["abundance"] = df["cov_sim"] / df["sum_cov"]

    if "reads" in entry_df.columns:
        df["reads"] = df["reads"] / df["count"]
        df["bases"] = df["reads"] * (snakemake.params.read_len * p)
        df["sum_bases"] = df.groupby("samplename")["bases"].transform("sum")
        df["proportion"] = df["bases"] / df["sum_bases"]
        df["cov_sim"] = df["bases"] / df["genome_size"]
        df["sum_cov"] = df.groupby("samplename")["cov_sim"].transform("sum")
        df["abundance"] = df["cov_sim"] / df["sum_cov"]

    if "bases" in entry_df.columns:
        df["bases"] = df["bases"] / df["count"]
        df["reads"] = df["bases"] / (snakemake.params.read_len * p)
        df["sum_bases"] = df.groupby("samplename")["bases"].transform("sum")
        df["proportion"] = df["bases"] / df["sum_bases"]
        df["cov_sim"] = df["bases"] / df["genome_size"]
        df["sum_cov"] = df.groupby("samplename")["cov_sim"].transform("sum")
        df["abundance"] = df["cov_sim"] / df["sum_cov"]

    elif "cov_sim" in entry_df.columns:
        df["cov_sim"] = df["cov_sim"] / df["count"]
        df["sum_cov"] = df.groupby("samplename")["cov_sim"].transform("sum")
        df["abundance"] = df["cov_sim"] / df["sum_cov"]
        df["bases"] = df["cov_sim"] * df["genome_size"]
        df["sum_bases"] = df.groupby("samplename")["bases"].transform("sum")
        df["reads"] = df["bases"] / (snakemake.params.read_len * p)
        df["proportion"] = df["bases"] / df["sum_bases"]


df["fasta"] = [str(fa).split("/")[-1].split(".fna.gz")[0] for fa in df["path"]]
df["seed"] = random.sample(range(1, 1000000), len(df))

cols = [
    "samplename",
    "entry",
    "fasta",
    "path",
    "taxid",
    "genome_size",
    "contig_count",
    "reads",
    "bases",
    "proportion",
    "cov_sim",
    "abundance",
    "seed",
]
df = df.astype({"seed": int})
df[cols].to_csv(snakemake.log[0], sep="\t", index=None)  # type: ignore
df[cols].set_index(["samplename", "fasta"]).to_csv(snakemake.output[0], sep="\t")  # type: ignore
