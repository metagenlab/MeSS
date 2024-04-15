from Bio import SeqIO
import pandas as pd
import os
from itertools import chain


def split_fasta(fa, outdir):
    record_ids = []
    name = os.path.basename(fa).split(".fasta")[0]
    fasta_name = os.path.join(outdir, name)

    for record in SeqIO.parse(fa, "fasta"):
        SeqIO.write(record, fasta_name + "_" + record.id + ".fna", "fasta")
        record_ids.append({"contig": record.id, "fasta": name})
    return record_ids


cov_df = pd.read_csv(snakemake.input.cov, sep="\t", dtype={"tax_id": int, "seed": int})

os.mkdir(snakemake.output.dir)
id2fa = []
for fa in snakemake.input.fa:
    id2fa.append(split_fasta(fa, snakemake.output.dir))
id2fa = list(chain.from_iterable(id2fa))
contig_df = pd.DataFrame.from_records(id2fa)
df = pd.merge(contig_df, cov_df, how="left", on="fasta")
df[["samplename", "fasta", "contig", "tax_id", "seed", "cov_sim"]].to_csv(
    snakemake.output.tsv, sep="\t", index=None
)
