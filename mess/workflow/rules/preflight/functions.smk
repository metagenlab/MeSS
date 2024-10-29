import os
import pandas as pd
import glob
from itertools import chain
from itertools import product
from Bio import SeqIO
import random


wildcard_constraints:
    sample="[^/]+",
    contig="[^/]+",
    fasta="[^/]+",


def list_reads(wildcards):
    fastqs = "{sample}.fq.gz"
    args = {"sample": SAMPLES}

    if PAIRED:
        fastqs = "{sample}_R{p}.fq.gz"
        args.update({"p": PAIRS})
    reads = collect(os.path.join(dir.out.fastq, fastqs), **args)

    if BAM:
        bams = collect(
            os.path.join(dir.out.bam, "{sample}.{bam}"),
            sample=SAMPLES,
            bam=["bam", "bam.bai"],
        )
        tax = collect(
            os.path.join(dir.out.tax, "{sample}_{abundance}.txt"),
            sample=SAMPLES,
            abundance=["seq", "tax"],
        )
        reads = reads + bams + tax

    return reads


def parse_samples(indir, replicates):
    if os.path.isfile(indir):
        files = [indir]
    else:
        files = glob.glob(f"{indir}/*.tsv")

    samples = []
    for file in files:
        try:
            samples.append(
                list(pd.read_csv(file, sep="\t")["sample"].drop_duplicates())
            )
        except KeyError:
            samples.append(os.path.basename(file).split(".")[0])
    if any(isinstance(sample, list) for sample in samples):
        samples = list(chain.from_iterable(samples))
    if len(replicates) < 2:
        return samples
    else:
        products = list(product(samples, replicates))
        return [prod[0] + "-" + str(prod[1]) for prod in products]


fasta_cache = {}


def get_fasta_table(wildcards):
    fa_table = checkpoints.calculate_genome_coverages.get(**wildcards).output[0]
    if fa_table not in fasta_cache:
        fa_df = pd.read_csv(fa_table, sep="\t", index_col="fasta")
        fasta_cache[fa_table] = fa_df
    fa_df = fasta_cache[fa_table]
    return fa_df


def fasta_input(wildcards):
    df = get_fasta_table(wildcards)
    try:
        return df.loc[wildcards.fasta]["path"].drop_duplicates()
    except AttributeError:
        return df.loc[wildcards.fasta]["path"]


def list_fastas(wildcards):
    df = get_fasta_table(wildcards)
    return expand(
        os.path.join(dir.out.processing, "{fasta}.fasta"), fasta=list(set(df.index))
    )


table_cache = {}


def get_cov_table(wildcards, key, idx_col):
    cov_table = checkpoints.split_contigs.get(**wildcards).output[0]
    if cov_table not in table_cache:
        cov_df = pd.read_csv(cov_table, sep="\t", index_col=idx_col).sort_index()
        table_cache[key] = cov_df
    cov_df = table_cache[key]
    return cov_df


def get_value(value, wildcards):
    vals = (
        f"{wildcards.sample}",
        f"{wildcards.fasta}",
        f"{wildcards.contig}",
    )
    idx_col = ["samplename", "fasta", "contig"]

    if CIRCULAR:
        idx_col += ["n"]
        vals += (int(wildcards.n),)
    val_df = get_cov_table(wildcards, "values", idx_col)
    return val_df.loc[vals, value]


def get_asm_summary(wildcards):
    try:
        table = checkpoints.download_assemblies.get(**wildcards).output[0]

    except AttributeError:
        if FASTA and not ASM_SUMMARY:
            table = os.path.join(dir.out.processing, "seqkit_stats.tsv")
        else:
            table = ASM_SUMMARY
    return table


def is_circular():
    if os.path.isfile(INPUT):
        files = [INPUT]
    else:
        files = glob.glob(f"{INPUT}/*.tsv")
    df = pd.concat([pd.read_csv(file, sep="\t") for file in files])
    if ROTATE > 1 or "rotate" in df.columns:
        return True
    else:
        return False


def aggregate(wildcards, outdir, ext):
    df = get_cov_table(wildcards, "aggregate", ["samplename"])
    files = []
    for row in df.loc[[wildcards.sample]].itertuples():
        prefix = f"{row.contig}"
        if CIRCULAR:
            prefix += f"_{row.n}"
        if PAIRED and ext != "bam":
            prefix += f"{wildcards.p}"
        files.append(
            os.path.join(outdir, wildcards.sample, row.fasta, f"{prefix}.{ext}")
        )
    return files


def get_header(fa):
    seqid = SeqIO.read(fa, "fasta").id
    return seqid.replace(">", "")


def seqkit_replace(wildcards):
    if PAIRED:
        return "{sample}_{nr}/{p}".format(
            sample=wildcards.sample, nr="{nr}", p=wildcards.p
        )
    else:
        return "{sample}_{nr}".format(
            sample=wildcards.sample,
            nr="{nr}",
        )
