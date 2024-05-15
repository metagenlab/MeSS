import pandas as pd
import numpy as np


# Set seed for replicates
np.random.seed(snakemake.params.seed)
# Duplicate each sample by the number of replicates

df = pd.read_csv(snakemake.input[0], sep="\t")
samples_df = df.reindex(df.index.repeat(len(snakemake.params.replicates)))
# add replicate column
if len(snakemake.params.replicates) < 2:
    samples_df.rename(columns={"sample": "samplename"}, inplace=True)
else:
    samples_df["rep"] = snakemake.params.replicates * len(samples_df)
    samples_df["samplename"] = [
        sample + "-" + str(rep)
        for sample, rep in zip(samples_df["sample"], samples_df["rep"])
    ]
    samples_df.drop(["sample", "rep"], axis=1, inplace=True)
accepted_cols = ["tax_abundance", "seq_abundance", "reads", "bases", "cov_sim"]

for col in accepted_cols:
    try:
        samples_df[f"{col}"] = [
            np.random.normal(val, snakemake.params.rep_sd)
            for val in samples_df[f"{col}"]
        ]
    except KeyError:
        continue

samples_df.reset_index(inplace=True, drop=True)
samples_df.to_csv(snakemake.output[0], sep="\t", index=None)
