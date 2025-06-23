from Bio import SeqIO
import pandas as pd
import os
from itertools import chain
import random
from snakemake.script import snakemake


def get_random_start(seed, contig_length, n):
    random.seed(seed)
    return [random.randint(0, contig_length - 1) for _ in range(n)]


def split_fasta(fa, outdir, suffix):
    record_ids = []
    name = os.path.basename(fa).split(suffix)[0]
    fasta_name = os.path.join(outdir, name)
    for record in SeqIO.parse(fa, "fasta"):
        SeqIO.write(record, fasta_name + "_" + record.id + ".fna", "fasta")
        record_ids.append(
            {"contig": record.id, "fasta": name, "contig_length": len(record.seq)}
        )
    return record_ids


cov_df = (
    pd.read_csv(snakemake.input.cov, sep="\t").dropna().convert_dtypes()
)  # dropna for dropping empty amplicons fasta


os.mkdir(snakemake.output.dir)
id2fa = []
suffix = ".fasta"
if snakemake.params.amplicons:
    suffix = ".amplicons.fasta"
for fa in snakemake.input.fa:
    id2fa.append(split_fasta(fa, snakemake.output.dir, suffix))
id2fa = list(chain.from_iterable(id2fa))
contig_df = pd.DataFrame.from_records(id2fa)

df = pd.merge(contig_df, cov_df, how="left", on="fasta")

cols = [
    "samplename",
    "fasta",
    "contig",
    "contig_length",
    "cov_sim",
]

if snakemake.params.circular:
    cols += ["n", "random_start", "rotate"]
    if "rotate" not in df.columns:
        df.loc[:, "rotate"] = [snakemake.params.rotate] * len(df)
    df["random_start"] = df.apply(
        lambda row: get_random_start(row["seed"], row["contig_length"], row["rotate"]),
        axis=1,
    )
    df_expanded = df.explode("random_start").reset_index(drop=True)
    df_expanded["n"] = df_expanded.groupby(["samplename", "contig"]).cumcount() + 1
    df = df_expanded
    df["cov_sim"] = df["cov_sim"] / df["rotate"]


cols += ["reads", "bases", "seq_abundance", "tax_abundance", "seed"]

if "tax_id" in df.columns:
    cols.append("tax_id")

df[cols].to_csv(snakemake.output.tsv, sep="\t", index=False)
