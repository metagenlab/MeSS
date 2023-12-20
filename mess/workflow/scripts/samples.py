#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
import glob

if os.path.isfile(snakemake.params.input):
    files = snakemake.params.input
else:
    files = glob.glob(f"{snakemake.params.input}/*.tsv")
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

# Set seed for replicates
np.random.seed(snakemake.params.seed)
# Duplicate each sample by the number of replicates

samples_df = dfs.reindex(dfs.index.repeat(len(snakemake.params.replicates)))
# add replicate column
samples_df["rep"] = snakemake.params.replicates * len(samples_df)
samples_df["samplename"] = [
    sample + "-" + str(rep)
    for sample, rep in zip(samples_df["sample"], samples_df["rep"])
]

samples_df.drop(["sample", "rep"], axis=1, inplace=True)
accepted_cols = ["proportion", "reads", "bases", "cov_sim"]

for col in accepted_cols:
    try:
        samples_df[f"{col}"] = [
            np.random.normal(val, snakemake.params.rep_sd) for val in samples_df[f"{col}"]
        ]
    except KeyError:
        continue

samples_df.reset_index(inplace=True, drop=True)
samples_df = samples_df.reindex(samples_df.index.repeat(len(snakemake.params.chunks)))
samples_df["chunk"] = snakemake.params.chunks * len(samples_df)
samples_df.to_csv(snakemake.output[0], sep="\t", index=None)
