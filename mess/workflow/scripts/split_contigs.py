from Bio import SeqIO
import pandas as pd
import os
from itertools import chain
import random


def get_random_start(seed, contig_length, n):
    random.seed(seed)
    return [random.randint(1, contig_length) for _ in range(n)]


def split_fasta(fa, outdir):
    record_ids = []
    name = os.path.basename(fa).split(".fasta")[0]
    fasta_name = os.path.join(outdir, name)

    for record in SeqIO.parse(fa, "fasta"):
        SeqIO.write(record, fasta_name + "_" + record.id + ".fna", "fasta")
        record_ids.append(
            {"contig": record.id, "fasta": name, "contig_length": len(record.seq)}
        )
    return record_ids


cov_df = pd.read_csv(snakemake.input.cov, sep="\t", dtype={"tax_id": int, "seed": int})

os.mkdir(snakemake.output.dir)
id2fa = []
for fa in snakemake.input.fa:
    id2fa.append(split_fasta(fa, snakemake.output.dir))
id2fa = list(chain.from_iterable(id2fa))
contig_df = pd.DataFrame.from_records(id2fa)
df = pd.merge(contig_df, cov_df, how="left", on="fasta")

cols = ["samplename", "fasta", "contig", "contig_length", "tax_id", "seed", "cov_sim"]

if snakemake.params.rotate > 1:
    cols += ["n", "contig_start"]
    df["contig_start"] = df.apply(
        lambda row: get_random_start(
            row["seed"], row["contig_length"], snakemake.params.rotate
        ),
        axis=1,
    )
    df_expanded = df.explode("contig_start").reset_index(drop=True)
    df_expanded["n"] = df_expanded.groupby(["samplename", "contig"]).cumcount() + 1
    df = df_expanded
    df["cov_sim"] = df["cov_sim"] / snakemake.params.rotate

df[cols].to_csv(snakemake.output.tsv, sep="\t", index=False)
