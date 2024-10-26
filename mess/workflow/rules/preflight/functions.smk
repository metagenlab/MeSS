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


def list_reads(wildcards):
    if PAIRED:
        reads = expand(
            os.path.join(dir.out.fastq, "{sample}_R{p}.fq.gz"),
            sample=SAMPLES,
            p=PAIRS,
        )
    else:
        reads = expand(
            os.path.join(dir.out.fastq, "{sample}.fq.gz"),
            sample=SAMPLES,
        )

    if BAM:
        bams = expand(
            os.path.join(dir.out.bam, "{sample}.{bam}"),
            sample=SAMPLES,
            bam=["bam", "bam.bai"],
        )
        tax = expand(
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


def fasta_input(wildcards):
    table = checkpoints.calculate_genome_coverages.get(**wildcards).output[0]

    df = pd.read_csv(table, sep="\t", index_col="fasta")
    try:
        return df.loc[wildcards.fasta]["path"].drop_duplicates()
    except AttributeError:
        return df.loc[wildcards.fasta]["path"]
    # some samples use the same genome path, drop duplicates to avoid duplicate paths when processing fasta


def list_fastas(wildcards):
    table = checkpoints.calculate_genome_coverages.get(**wildcards).output[0]
    if table not in fasta_cache:
        df = pd.read_csv(table, sep="\t")
        fasta_cache[table] = df
    df = fasta_cache[table]
    fastas = list(set(df["fasta"]))
    return expand(os.path.join(dir.out.processing, "{fasta}.fasta"), fasta=fastas)


table_cache = {}


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

    table = checkpoints.split_contigs.get(**wildcards).output[0]
    if table not in table_cache:
        df = pd.read_csv(
            table,
            sep="\t",
            index_col=idx_col,
        ).sort_index()
        table_cache[table] = df
    df = table_cache[table]
    return df.loc[vals, value]


def get_asm_summary(wildcards):
    try:
        table = checkpoints.download_assemblies.get(**wildcards).output[0]

    except AttributeError:
        if FASTA and not ASM_SUMMARY:
            table = os.path.join(dir.out.processing, "seqkit_stats.tsv")
        else:
            table = ASM_SUMMARY
    return table


def get_cov_table(wildcards):
    return checkpoints.split_contigs.get(**wildcards).output[0]


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
    table = checkpoints.split_contigs.get(**wildcards).output[0]
    df = pd.read_csv(
        table,
        sep="\t",
        index_col=["samplename", "fasta"],
    ).sort_index()
    fastas = list(set(df.loc[wildcards.sample].index))
    contigs = list(
        chain(*[list(df.loc[(wildcards.sample, fasta), "contig"]) for fasta in fastas])
    )

    collect_args = {
        "sample": wildcards.sample,
        "fasta": fastas,
        "contig": contigs,
        "ext": ext,
    }
    path = os.path.join(outdir, "{sample}", "{fasta}", "{contig}.{ext}")
    if CIRCULAR:
        path = os.path.join(outdir, "{sample}", "{fasta}", "{contig}_{n}.{ext}")
        rotates = list(
            chain(*[list(df.loc[(wildcards.sample, fasta), "n"]) for fasta in fastas])
        )
        collect_args.update(
            {
                "n": rotates,
            }
        )
    if PAIRED and ext != "bam":
        path = os.path.join(outdir, "{sample}", "{fasta}", "{contig}{p}.{ext}")
        collect_args.update(
            {
                "p": wildcards.p,
            }
        )
        if CIRCULAR:
            path = os.path.join(outdir, "{sample}", "{fasta}", "{contig}_{n}{p}.{ext}")
    return collect(path, **collect_args)


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
