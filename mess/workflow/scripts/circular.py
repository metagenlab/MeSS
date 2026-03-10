from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import random
from snakemake.script import snakemake
import os
import re


def strip_fasta_ext(filename):
    return re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", filename)


df = pd.read_csv(snakemake.input.fa2id, sep="\t", header=None)
df.columns = ["fasta", "seq_id", "length"]
pattern = "|".join(snakemake.params.circular.split(","))
matches = df[df["seq_id"].str.contains(pattern, regex=True)]
rotates = snakemake.params.rotates
random.seed(snakemake.params.seed)
for fasta in matches.fasta.drop_duplicates().to_list():
    basename = strip_fasta_ext(os.path.basename(fasta))
    for record in SeqIO.parse(
        os.path.join(snakemake.params.fadir, f"{basename}.fasta"), "fasta"
    ):
        if matches["seq_id"].str.contains(record.id).any():
            ori_seq = record.seq
            seq_id = record.id
            records = []
            for n in range(rotates):
                start = random.randint(0, len(ori_seq) - 1)
                end = start + len(ori_seq)
                records.append(
                    SeqRecord(
                        ori_seq[start:] + ori_seq[:start],
                        id=f"{seq_id}_circular_{n + 1}_{start}:{end}",
                        description="",
                    )
                )
            SeqIO.write(
                records,
                os.path.join(snakemake.params.fadir, f"{basename}_circular.fasta"),
                "fasta",
            )
        else:
            record.description = ""
            SeqIO.write(
                record,
                os.path.join(snakemake.params.fadir, f"{basename}.fasta"),
                "fasta",
            )
with open(snakemake.output[0], "w") as f:
    f.write("circular.done\n")
