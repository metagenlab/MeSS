from Bio import SeqIO
import pandas as pd
import os
from itertools import chain


def split_fasta(fa, outdir):
    record_ids = []
    name = os.path.basename(fa).split(".fasta")[0]
    os.mkdir(os.path.join(outdir, name))
    for record in SeqIO.parse(fa, "fasta"):
        SeqIO.write(record, os.path.join(outdir, name, record.id) + ".fa", "fasta")
        record_ids.append({"contig": record.id, "fasta": name})
    return record_ids


id2fa = []
for fa in snakemake.input:
    id2fa.append(split_fasta(fa, snakemake.params[0]))
id2fa = list(chain.from_iterable(id2fa))
pd.DataFrame.from_records(id2fa).to_csv(snakemake.output[0], sep="\t", index=None)
