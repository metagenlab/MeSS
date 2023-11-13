#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import logging
import os

"""
Description
"""
desc = "calculate reads, bases, coverage, proportions and abundances"
example = """Example usage: get_dist.py -i input.tsv -a assembly_summary.tsv -o abundances.json"""
argParser = argparse.ArgumentParser(description=desc, epilog=example)

"""
Arguments
"""
argParser.add_argument(
    "-i",
    dest="input",
    type=str,
    help="mess input table",
    required=True,
)
argParser.add_argument(
    "-a",
    dest="asm_summary",
    type=str,
    help="assembly summary table with taxonomy and sequence info",
    required=True,
)
argParser.add_argument(
    "-o", dest="output", type=str, help="output json", default="abundances.json"
)

argParser.add_argument(
    "-d",
    dest="dist",
    type=list,
    help="distribution of reads across genomes",
    choices=["lognormal", "even", None],
    default=None,
)

argParser.add_argument(
    "-m",
    "--mu",
    dest="mu",
    type=int,
    help="mu of lognormal distribution",
    default=1,
)
argParser.add_argument(
    "-s",
    "--sigma",
    dest="sigma",
    type=int,
    help="sigma of lognormal distribution",
    default=0,
)

argParser.add_argument(
    "--log", dest="log", type=str, help="log path", default="log.txt"
)
argParser.add_argument(
    "-n",
    "--total_bases",
    dest="total_bases",
    type=str,
    help="total sequencing bases",
    default="1G",
)

argParser.add_argument(
    "-l",
    "--read_length",
    dest="read_length",
    type=int,
    help="average read length",
    default=150,
)

argParser.add_argument(
    "-p",
    "--pairing",
    dest="pairing",
    type=bool,
    help="are reads paired or not",
    default=False,
)

argParser.add_argument(
    "-sd",
    "--sd_rep",
    dest="sd_rep",
    type=int,
    help="standard deviation of read count between replicates",
    default=0,
)
argParser.add_argument(
    "--seed",
    dest="seed",
    type=int,
    help="set seed for reproducibility",
    default=42,
)

args = argParser.parse_args()

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


def conv(s):
    try:
        # multiply with meter-prefix value
        return float(s[:-1]) * prefix[s[-1]]
    except KeyError:
        # no or unknown meter-prefix
        return float(s)


"""
Main
"""
# Set logging
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%d %b %Y %H:%M:%S",
    filename=args.log,
    level=logging.DEBUG,
)
# Set seed
np.random.seed(args.seed)

# Pairing
if args.pairing:
    p = 2
else:
    p = 1

# Calculate total base count

prefix = {"T": 10**12, "G": 10**9, "M": 10**6, "k": 10**3}
total_bases = conv(args.total_bases)


# Get table with assembly genomsizes and their taxonomy
entry_df = pd.read_csv(args.input, sep="\t")
asm_df = pd.read_csv(args.asm_summary, sep="\t")
df = entry_df.merge(asm_df)

# Count number of entries
df["count"] = df.groupby("entry")["entry"].transform("count")

# Calculate prportion with dist
if args.dist == "even":
    df = get_even_dist(df, ["taxid"])
    df["bases"] = df["proportion"] * total_bases
    df["reads"] = df["bases"] / (args.read_length * p)
    df["cov_sim"] = df["bases"] / df["genome_size"]
    df["abundance"] = df["cov_sim"] / df["cov_sim"].sum()

elif args.dist == "lognormal":
    df = get_lognormal_dist(df, args.mu, args.sigma)
    df["bases"] = df["proportion"] * total_bases
    df["reads"] = df["bases"] / (args.read_length * p)
    df["cov_sim"] = df["bases"] / df["genome_size"]
    df["abundance"] = df["cov_sim"] / df["cov_sim"].sum()

elif args.dist == None:
    if "proportion" in entry_df.columns:
        df["proportion"] = df["proportion"] / df["count"]
        df["bases"] = df["proportion"] * total_bases
        df["reads"] = df["bases"] / (args.read_length * p)
        df["cov_sim"] = df["bases"] / df["genome_size"]
        df["abundance"] = df["cov_sim"] / df["cov_sim"].sum()

    elif "reads" in entry_df.columns:
        df["reads"] = df["reads"] / df["count"]
        df["bases"] = df["reads"] * (args.read_length * p)
        total_bases = df["bases"].sum()
        df["proportion"] = df["bases"] / total_bases
        df["cov_sim"] = df["bases"] / df["genome_size"]
        df["abundance"] = df["cov_sim"] / df["cov_sim"].sum()

    elif "bases" in entry_df.columns:
        total_bases = df["bases"].sum()
        df["bases"] = df["bases"] / df["count"]
        df["reads"] = df["bases"] / (args.read_length * p)
        df["proportion"] = df["bases"] / total_bases
        df["cov_sim"] = df["bases"] / df["genome_size"]
        df["abundance"] = df["cov_sim"] / df["cov_sim"].sum()

    elif "cov_sim" in entry_df.columns:
        df["cov_sim"] = df["cov_sim"] / df["count"]
        df["bases"] = df["cov_sim"] * df["genome_size"]
        df["reads"] = df["bases"] / (args.read_length * p)
        df["proportion"] = df["bases"] / df["bases"].sum()
        df["abundance"] = df["cov_sim"] / df["cov_sim"].sum()

df["fasta"] = df["path"].apply(os.path.basename)
df["cov_sim"] = [np.random.normal(cov, args.sd_rep) for cov in df["cov_sim"]]
sel_df = df[["fasta", "proportion", "bases", "reads", "cov_sim", "abundance"]]

sel_df.to_csv(args.log, sep="\t", index=None)
sel_df.set_index("fasta").to_json(args.output, index=None)
