import pandas as pd
import numpy as np
from humanfriendly import parse_size
import random
import os
import glob
import re
from snakemake.script import snakemake

"""
Functions
"""


def get_lognormal_dist(df, mu, sigma):
    """
    function for simulating lognormal sequence distribution
    """
    dfs = []
    for sample in set(df["samplename"]):
        subdf = df[df["samplename"] == sample]
        subdf.loc[:, "lognormal"] = np.random.lognormal(
            mean=mu, sigma=sigma, size=len(subdf)
        )
        subdf.loc[:, "seq_abundance"] = subdf["lognormal"] / subdf["lognormal"].sum()
        dfs.append(subdf)
    return pd.concat(dfs)


def get_even_dist(table):
    """
    function that calculates even abundances across taxonomy
    """
    abundances = (
        table.groupby(["samplename"])["tax_id"]
        .value_counts(normalize=True)
        .reset_index()
    )
    counts = table.groupby("samplename")["tax_id"].value_counts().reset_index()
    return table.merge(abundances).merge(counts)


def strip_fasta_ext(filename):
    return re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", filename)


"""
Main
"""
# Set seed for distributions
np.random.seed(snakemake.params["seed"])

# Set seed for read simulators
random.seed(snakemake.params["seed"])

# Pairing
if snakemake.params["pairing"]:
    p = 2
else:
    p = 1

entry_df = pd.read_csv(snakemake.input["rep"], sep="\t")


if ("fasta" not in entry_df.columns) and ("path" in entry_df.columns):
    entry_df["fasta"] = [
        strip_fasta_ext(os.path.basename(path)) for path in entry_df["path"]
    ]


if snakemake.params["fa_dir"]:
    entry_df["fasta"] = entry_df["fasta"].apply(strip_fasta_ext)
    entry_df["path"] = [
        glob.glob(os.path.join(snakemake.params["fa_dir"], f"{fa}*"))[0]
        for fa in entry_df["fasta"]
    ]
if snakemake.params["fa_path"]:
    entry_df["fasta"] = [
        strip_fasta_ext(os.path.basename(path)) for path in entry_df["path"]
    ]

if isinstance(snakemake.input["asm"], list):
    summary_df = pd.read_csv(snakemake.input["asm"][0], sep="\t")
    summary_df["fasta"] = [
        strip_fasta_ext(os.path.basename(path)) for path in summary_df["path"]
    ]
    stats_df = pd.read_csv(snakemake.input["asm"][1], sep="\t")
    stats_df.rename(
        columns={
            "file": "fasta",
            "sum_len": "seq_len",
            "num_seqs": "seq_num",
        },
        inplace=True,
    )
    stats_df["fasta"] = stats_df["fasta"].apply(strip_fasta_ext)
    stats_df["fasta"] = stats_df["fasta"].str.replace(".amplicons", "")
    asm_df = pd.merge(stats_df, summary_df, on="fasta")


else:
    asm_df = pd.read_csv(snakemake.input["asm"], sep="\t")
    if "file" in asm_df.columns:
        asm_df.rename(
            columns={
                "file": "fasta",
                "sum_len": "seq_len",
                "num_seqs": "seq_num",
            },
            inplace=True,
        )
        asm_df["fasta"] = asm_df["fasta"].apply(strip_fasta_ext)

    if ("fasta" not in asm_df.columns) and ("path" in asm_df.columns):
        asm_df["fasta"] = [
            strip_fasta_ext(os.path.basename(path)) for path in asm_df["path"]
        ]
    if snakemake.params["amplicons"]:
        asm_df["fasta"] = asm_df["fasta"].str.replace(".amplicons", "")
    if (
        "total_sequence_length" in asm_df.columns
        and "number_of_contigs" in asm_df.columns
    ):
        asm_df = asm_df.rename(
            columns={"total_sequence_length": "seq_len", "number_of_contigs": "seq_num"}
        )


same_cols = list(np.intersect1d(entry_df.columns, asm_df.columns))
df = pd.merge(entry_df, asm_df, how="left", on=same_cols)


# Get total bases
bases = parse_size(snakemake.params["bases"])


if "tax_id" in df.columns:
    df["tax_id"] = df["tax_id"].astype(int)
# Calculate prportion with dist
if snakemake.params["dist"] == "even":
    df = get_even_dist(df)
    df["tax_abundance"] = df["proportion"] / df["count"]
    df["genome_bases"] = df["seq_len"] * df["tax_abundance"]
    df["sum_genome_bases"] = df.groupby("samplename")["genome_bases"].transform("sum")
    df["cov_obtained"] = bases / df["sum_genome_bases"]
    df["cov_sim"] = df["tax_abundance"] * df["cov_obtained"]
    df["sum_cov"] = df.groupby("samplename")["cov_sim"].transform("sum")
    df["bases"] = df["cov_sim"] * df["seq_len"]
    df["reads"] = df["bases"] / snakemake.params["read_len"]
    df["sum_bases"] = df.groupby("samplename")["bases"].transform("sum")
    df["seq_abundance"] = df["bases"] / df["sum_bases"]


elif snakemake.params["dist"] == "lognormal":
    df = get_lognormal_dist(
        df, mu=snakemake.params["mu"], sigma=snakemake.params["sigma"]
    )
    df["bases"] = df["seq_abundance"] * bases
    df["reads"] = df["bases"] / snakemake.params["read_len"]
    df["cov_sim"] = df["bases"] / df["seq_len"]
    df["sum_cov"] = df.groupby("samplename")["cov_sim"].transform("sum")
    df["tax_abundance"] = df["cov_sim"] / df["sum_cov"]
    df["sum_bases"] = df.groupby("samplename")["bases"].transform("sum")
    df["seq_abundance"] = df["bases"] / df["sum_bases"]
else:
    if "tax_abundance" in entry_df.columns:
        df["genome_bases"] = df["seq_len"] * df["tax_abundance"]
        df["sum_genome_bases"] = df.groupby("samplename")["genome_bases"].transform(
            "sum"
        )
        df["cov_obtained"] = bases / df["sum_genome_bases"]
        df["cov_sim"] = df["tax_abundance"] * df["cov_obtained"]
        df["sum_cov"] = df.groupby("samplename")["cov_sim"].transform("sum")
        df["bases"] = df["cov_sim"] * df["seq_len"]
        df["reads"] = df["bases"] / snakemake.params["read_len"]
        df["sum_bases"] = df.groupby("samplename")["bases"].transform("sum")
        df["seq_abundance"] = df["bases"] / df["sum_bases"]

    if "seq_abundance" in entry_df.columns:
        df["bases"] = df["seq_abundance"] * bases
        df["reads"] = df["bases"] / snakemake.params["read_len"]
        df["cov_sim"] = df["bases"] / df["seq_len"]
        df["sum_cov"] = df.groupby("samplename")["cov_sim"].transform("sum")
        df["tax_abundance"] = df["cov_sim"] / df["sum_cov"]

    if "reads" in entry_df.columns:
        df["bases"] = df["reads"] * snakemake.params["read_len"] * p
        df["sum_bases"] = df.groupby("samplename")["bases"].transform("sum")
        df["seq_abundance"] = df["bases"] / df["sum_bases"]
        df["cov_sim"] = df["bases"] / df["seq_len"]
        df["sum_cov"] = df.groupby("samplename")["cov_sim"].transform("sum")
        df["tax_abundance"] = df["cov_sim"] / df["sum_cov"]

    if "bases" in entry_df.columns:
        df["reads"] = df["bases"] / snakemake.params["read_len"]
        df["sum_bases"] = df.groupby("samplename")["bases"].transform("sum")
        df["seq_abundance"] = df["bases"] / df["sum_bases"]
        df["cov_sim"] = df["bases"] / df["seq_len"]
        df["sum_cov"] = df.groupby("samplename")["cov_sim"].transform("sum")
        df["tax_abundance"] = df["cov_sim"] / df["sum_cov"]

    elif "cov_sim" in entry_df.columns:
        df["sum_cov"] = df.groupby("samplename")["cov_sim"].transform("sum")
        df["tax_abundance"] = df["cov_sim"] / df["sum_cov"]
        df["bases"] = df["cov_sim"] * df["seq_len"]
        df["sum_bases"] = df.groupby("samplename")["bases"].transform("sum")
        df["reads"] = df["bases"] / snakemake.params["read_len"]
        df["seq_abundance"] = df["bases"] / df["sum_bases"]


df["seed"] = random.sample(range(1, 1000000), len(df))
cols = [
    "samplename",
    "fasta",
    "seq_len",
    "seq_num",
    "reads",
    "bases",
    "seq_abundance",
    "cov_sim",
    "tax_abundance",
    "seed",
]
if "rotate" in entry_df.columns:
    cols.append("rotate")
if "tax_id" in df.columns:
    cols.append("tax_id")

# replace values with 0 for empty amplicon fastas

if (df["seq_len"] == 0).any():
    df.loc[
        df["seq_len"] == 0,
        [
            "seq_num",
            "reads",
            "bases",
            "cov_sim",
            "tax_abundance",
            "seq_abundance",
            "seed",
        ],
    ] = 0
    df = df[cols].replace(0, np.nan)
df[cols].sort_values(
    ["samplename", "fasta", "cov_sim"], ascending=[True, True, False]
).convert_dtypes().to_csv(snakemake.output[0], sep="\t", index=False)
