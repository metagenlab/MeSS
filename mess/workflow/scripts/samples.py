#!/usr/bin/env python3

import os
import glob
import pandas as pd

if os.path.isfile(snakemake.params[0]):
    files = snakemake.params[0]
else:
    files = glob.glob(f"{snakemake.params[0]}/*.tsv")

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
dfs.to_csv(snakemake.output[0], sep="\t", index=None)
