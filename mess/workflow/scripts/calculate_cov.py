import pandas as pd
import numpy as np
from humanfriendly import parse_size
import random
import os
import glob

"""
Functions
"""


def get_lognormal_dist(df, mu, sigma):
    """
    function for simulating lognormal sequence distribution
    """
    df["lognormal"] = np.random.lognormal(mean=mu, sigma=sigma, size=len(df))
    df["tax_abundance"] = df["lognormal"] / df["lognormal"].sum()
    return df


def get_even_dist(df, cols):
    """
    function that calculates even abundances across taxonomy
    """
    abundances = (
        df[cols]
        .value_counts(normalize=True)
        .reset_index()
        .rename(columns={"seq_abundance": "tax_abundance"}, inplace=True)
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

if snakemake.params.fa:
    entry_df = pd.read_csv(snakemake.input.df, sep="\t")
    entry_df["path"] = [
        glob.glob(os.path.join(snakemake.params.fa, f"{fa}*"))[0]
        for fa in entry_df["fasta"]
    ]
    asm_df = pd.read_csv(snakemake.input.asm, sep="\t")
    asm_df.rename(
        columns={
            "file": "path",
            "sum_len": "total_sequence_length",
            "num_seqs": "number_of_contigs",
        },
        inplace=True,
    )

else:
    asm_df = pd.read_csv(snakemake.input.asm, sep="\t")
    entry_df = pd.read_csv(snakemake.input.df, sep="\t")


same_cols = list(np.intersect1d(entry_df.columns, asm_df.columns))
df = pd.merge(entry_df, asm_df, how="left", on=same_cols)

# Get total bases
bases = parse_size(snakemake.params.bases)

# Calculate prportion with dist
if snakemake.params.dist == "even":
    df = get_even_dist(df, ["tax_id"])
    df["sum_seq_length"] = df.groupby("samplename")["total_sequence_length"].transform(
        "sum"
    )
    df["sum_cov"] = bases / df["sum_seq_length"]
    df["cov_sim"] = df["sum_cov"] * df["tax_abundance"]


elif snakemake.params.dist == "lognormal":
    df = get_lognormal_dist(df, snakemake.params.mu, snakemake.params.sigma)
    df["bases"] = df["seq_abundance"] * bases
    df["reads"] = df["bases"] / (snakemake.params.read_len * p)
    df["cov_sim"] = df["bases"] / df["total_sequence_length"]
    df["sum_cov"] = df.groupby("samplename")["cov_sim"].transform("sum")
    df["tax_abundance"] = df["cov_sim"] / df["sum_cov"]
    df["sum_bases"] = df.groupby("samplename")["bases"].transform("sum")
    df["reads"] = df["bases"] / (snakemake.params.read_len * p)
    df["seq_abundance"] = df["bases"] / df["sum_bases"]
else:
    if "tax_abundance" in entry_df.columns:
        df["sum_seq_length"] = df.groupby("samplename")[
            "total_sequence_length"
        ].transform("sum")
        df["sum_cov"] = bases / df["sum_seq_length"]
        df["cov_sim"] = df["sum_cov"] * df["tax_abundance"]
        df["bases"] = df["cov_sim"] * df["total_sequence_length"]
        df["sum_bases"] = df.groupby("samplename")["bases"].transform("sum")
        df["seq_abundance"] = df["bases"] / df["sum_bases"]
        df["reads"] = df["bases"] / (snakemake.params.read_len * p)

    if "seq_abundance" in entry_df.columns:
        df["bases"] = df["seq_abundance"] * bases
        df["reads"] = df["bases"] / (snakemake.params.read_len * p)
        df["cov_sim"] = df["bases"] / df["total_sequence_length"]
        df["sum_cov"] = df.groupby("samplename")["cov_sim"].transform("sum")
        df["tax_abundance"] = df["cov_sim"] / df["sum_cov"]

    if "reads" in entry_df.columns:
        df["bases"] = df["reads"] * (snakemake.params.read_len * p)
        df["sum_bases"] = df.groupby("samplename")["bases"].transform("sum")
        df["seq_abundance"] = df["bases"] / df["sum_bases"]
        df["cov_sim"] = df["bases"] / df["total_sequence_length"]
        df["sum_cov"] = df.groupby("samplename")["cov_sim"].transform("sum")
        df["tax_abundance"] = df["cov_sim"] / df["sum_cov"]

    if "bases" in entry_df.columns:
        df["reads"] = df["bases"] / (snakemake.params.read_len * p)
        df["sum_bases"] = df.groupby("samplename")["bases"].transform("sum")
        df["seq_abundance"] = df["bases"] / df["sum_bases"]
        df["cov_sim"] = df["bases"] / df["total_sequence_length"]
        df["sum_cov"] = df.groupby("samplename")["cov_sim"].transform("sum")
        df["tax_abundance"] = df["cov_sim"] / df["sum_cov"]

    elif "cov_sim" in entry_df.columns:
        df["sum_cov"] = df.groupby("samplename")["cov_sim"].transform("sum")
        df["tax_abundance"] = df["cov_sim"] / df["sum_cov"]
        df["bases"] = df["cov_sim"] * df["total_sequence_length"]
        df["sum_bases"] = df.groupby("samplename")["bases"].transform("sum")
        df["reads"] = df["bases"] / (snakemake.params.read_len * p)
        df["seq_abundance"] = df["bases"] / df["sum_bases"]

if "fasta" not in df.columns:
    df["fasta"] = df["accession"]

df["seed"] = random.sample(range(1, 1000000), len(df))

cols = [
    "samplename",
    "fasta",
    "path",
    "tax_id",
    "total_sequence_length",
    "number_of_contigs",
    "reads",
    "bases",
    "seq_abundance",
    "cov_sim",
    "tax_abundance",
    "seed",
]
df = df.astype(
    {"seed": int, "tax_id": int, "total_sequence_length": int, "number_of_contigs": int}
)
df[cols].set_index(["samplename", "fasta"]).to_csv(snakemake.output[0], sep="\t")  # type: ignore
