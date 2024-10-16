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

    if ROTATE > 1:
        idx_col += ["n"]
        vals += (int(wildcards.n),)

    table = checkpoints.split_contigs.get(**wildcards).output[0]
    if table not in table_cache:
        df = pd.read_csv(
            table,
            sep="\t",
            index_col=idx_col,
        )
        table_cache[table] = df
    df = table_cache[table]
    return df.loc[vals, value]


fa_table_cache = {}
def get_random_start(wildcards):
    table = checkpoints.split_contigs.get(**wildcards).output[0]
    if table not in fa_table_cache:
        df = pd.read_csv(
            table,
            sep="\t",
            index_col=["fasta","contig","n"]
        )
        fa_table_cache[table] = df
    vals=(f"{wildcards.fasta}",f"{wildcards.contig}",int(wildcards.n))
    df = fa_table_cache[table]
    return df.loc[vals, "random_start"]





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


def aggregate(wildcards, outdir, level, ext):
    table = checkpoints.split_contigs.get(**wildcards).output[0]
    if level == "contig":
        df = pd.read_csv(table, sep="\t", index_col="fasta")
        contigs = df.loc[wildcards.fasta]["contig"]
        if isinstance(contigs, str):
            contigs = [contigs]
        else:
            contigs = list(set(contigs))
        if PAIRED and ext != "bam":
            if ROTATE > 1:
                return expand(
                    os.path.join(
                        outdir, "{sample}", "{fasta}", "{contig}_{n}{p}.{ext}"
                    ),
                    sample=wildcards.sample,
                    fasta=wildcards.fasta,
                    n=list(range(1, ROTATE + 1)),
                    p=wildcards.p,
                    contig=contigs,
                    ext=ext,
                )
            else:
                return expand(
                    os.path.join(outdir, "{sample}", "{fasta}", "{contig}{p}.{ext}"),
                    sample=wildcards.sample,
                    fasta=wildcards.fasta,
                    p=wildcards.p,
                    contig=contigs,
                    ext=ext,
                )

        else:
            if ROTATE > 1:
                return expand(
                    os.path.join(outdir, "{sample}", "{fasta}", "{contig}_{n}.{ext}"),
                    sample=wildcards.sample,
                    fasta=wildcards.fasta,
                    contig=contigs,
                    n=list(range(1, ROTATE + 1)),
                    ext=ext,
                )
            else:
                return expand(
                    os.path.join(outdir, "{sample}", "{fasta}", "{contig}.{ext}"),
                    sample=wildcards.sample,
                    fasta=wildcards.fasta,
                    contig=contigs,
                    ext=ext,
                )
    if level == "fasta":
        df = pd.read_csv(table, sep="\t", index_col=["samplename", "fasta"])
        fastas = list(set(df.loc[wildcards.sample].index))
        if PAIRED and ext != "bam":
            return expand(
                os.path.join(outdir, "{sample}", "{fasta}{p}.{ext}"),
                sample=wildcards.sample,
                fasta=fastas,
                p=wildcards.p,
                ext=ext,
            )
        else:
            return expand(
                os.path.join(outdir, "{sample}", "{fasta}.{ext}"),
                sample=wildcards.sample,
                fasta=fastas,
                ext=ext,
            )


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
