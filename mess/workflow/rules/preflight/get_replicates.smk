# samples and replicates tables

import pandas as pd
import os
import numpy as np
from itertools import chain
import glob


SEED = config.args.seed
# Set seed for replicates
np.random.seed(SEED)
# Read input path

# Duplicate each sample by the number of replicates
REPLICATES = list(range(1, config.args.replicates + 1))
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

REP_SD = config.args.rep_sd
for col in accepted_cols:
    try:
        SAMPLES_DF[f"{col}"] = [
            np.random.normal(val, REP_SD) for val in SAMPLES_DF[f"{col}"]
        ]
    except KeyError:
        continue
SAMPLES_DF.reset_index(inplace=True, drop=True)
CHUNKS = [f"{x+1:03}" for x in range(config.args.chunks)]
SAMPLES_DF = SAMPLES_DF.reindex(SAMPLES_DF.index.repeat(len(CHUNKS)))
SAMPLES_DF["chunk"] = CHUNKS * len(CHUNKS)


rule samples_table:
    input:
        os.path.join(dir.out.base, "cov.tsv"),


checkpoint calculate_coverage:
    input:
        asm=lambda wildcards: get_assembly_summary_path(wildcards),
    output:
        temp(
            os.path.join(dir.out.base, "cov.tsv"),
        ),
    params:
        dist=config.args.dist,
        mu=config.args.mu,
        sigma=config.args.sigma,
        total_bases=config.args.total_bases,
        read_len=config.args.read_len,
        pairing=config.args.paired,
        rep_sd=config.args.rep_sd,
        seed=config.args.seed,
        chunks=CHUNKS,
    log:
        os.path.join(dir.out.logs, "tables/cov.tsv"),
    script:
        "calculate_cov.py"
