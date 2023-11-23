#!/usr/bin/env python3

import random
import pandas as pd
import os
import numpy as np
from itertools import chain
import glob


def list_files(indir, extensions):
    extensions = extensions.split(",")
    files = [glob.glob(f"{indir}/*.{e}") for e in extensions]
    return list(chain.from_iterable(files))


# Set seed for replicates
np.random.seed(snakemake.params.seed)
# Read input path
if os.path.isfile(snakemake.input[0]):
    files = snakemake.input[0]
else:
    files = list_files(snakemake.input[0], "tsv")

# concat tables into one and add sample columns
dfs = []
for file in files:
    df = pd.read_csv(file, sep="\t")
    dfs.append(df)
    try:
        samples = list(set(df["sample"]))
    except KeyError:
        sample = os.path.basename(file).split(".")[0]
        df["sample"] = [sample] * len(df)

dfs = pd.concat(dfs).reset_index(drop=True).sort_values("sample")
# Duplicate each sample by the number of replicates
replicates = list(range(1, snakemake.params.rep + 1))
rep_df = dfs.reindex(dfs.index.repeat(len(replicates)))
# add replicate column
rep_df["rep"] = replicates * len(dfs)
rep_df["samplename"] = [
    sample + "-" + str(rep) for sample, rep in zip(rep_df["sample"], rep_df["rep"])
]
rep_df.drop(["sample", "rep"], axis=1, inplace=True)
accepted_cols = ["proportion", "reads", "bases", "cov_sim"]
for col in accepted_cols:
    try:
        rep_df[f"{col}"] = [
            np.random.normal(val, snakemake.params.rep_sd) for val in rep_df[f"{col}"]
        ]
    except KeyError:
        continue

rep_df.to_csv(snakemake.output[0], sep="\t", index=None)
