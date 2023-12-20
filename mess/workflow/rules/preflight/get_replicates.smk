# Get replicates, calculate coverage
import os
import numpy as np

# Set seed for replicates
np.random.seed(SEED)
# Duplicate each sample by the number of replicates
SAMPLES_DF = DFS.reindex(DFS.index.repeat(len(REPLICATES)))
# add replicate column
SAMPLES_DF["rep"] = REPLICATES * len(dfs)
SAMPLES_DF["samplename"] = [
    sample + "-" + str(rep)
    for sample, rep in zip(SAMPLES_DF["sample"], SAMPLES_DF["rep"])
]

SAMPLES_DF.drop(["sample", "rep"], axis=1, inplace=True)
SAMPLES = list(set(SAMPLES_DF["samplename"]))
accepted_cols = ["proportion", "reads", "bases", "cov_sim"]

for col in accepted_cols:
    try:
        SAMPLES_DF[f"{col}"] = [
            np.random.normal(val, REP_SD) for val in SAMPLES_DF[f"{col}"]
        ]
    except KeyError:
        continue

SAMPLES_DF.reset_index(inplace=True, drop=True)
SAMPLES_DF = SAMPLES_DF.reindex(SAMPLES_DF.index.repeat(len(CHUNKS)))
SAMPLES_DF["chunk"] = CHUNKS * len(CHUNKS)
