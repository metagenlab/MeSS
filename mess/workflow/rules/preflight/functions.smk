import os
import pandas as pd
import glob
import random
from itertools import chain
from itertools import product
from Bio import SeqIO


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
        tax = expand(os.path.join(dir.out.tax, "{sample}_profile.txt"), sample=SAMPLES)
        reads = reads + bams + tax

    return reads


def list_files(indir, extensions):
    extensions = extensions.split(",")
    files = [glob.glob(f"{indir}/*.{e}") for e in extensions]
    return list(chain.from_iterable(files))


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


def fasta_input(wildcards):
    table = checkpoints.calculate_coverage.get(**wildcards).output[0]
    df = pd.read_csv(table, sep="\t")
    basedir = os.path.dirname(df["path"][0])
    if COMPRESSED:
        return os.path.join(basedir, "{fasta}.fna.gz")
    elif COMPRESSED == False:
        return os.path.join(basedir, "{fasta}.fna")
    if SKIP_FA_PROC:
        return os.path.join(basedir, "{fasta}", "{contig}.fna")


def list_fastas(wildcards):
    table = checkpoints.calculate_coverage.get(**wildcards).output[0]
    df = pd.read_csv(table, sep="\t")
    fastas = list(set(df["fasta"]))
    return expand(os.path.join(dir.out.fasta, "{fasta}.fasta"), fasta=fastas)


table_cache = {}


def get_value(table, wildcards, value):
    if table not in table_cache:
        df = pd.read_csv(
            table,
            sep="\t",
            index_col=["samplename", "fasta"],
        )
        table_cache[table] = df
    df = table_cache[table]
    val = df.loc[wildcards.sample].loc[wildcards.fasta][value]
    return val


def get_assembly_summary_path(wildcards):
    try:
        table = os.path.abspath(
            checkpoints.download_assemblies.get(**wildcards).output[0]
        )
    except AttributeError:
        table = os.path.abspath(ASM_SUMMARY)
    return table


def aggregate(wildcards, outdir, level, ext):
    if level == "contig":
        table = checkpoints.split_contigs.get(**wildcards).output[0]
        df = pd.read_csv(table, sep="\t", index_col="fasta")
        contigs = df.loc[wildcards.fasta]["contig"]
        if isinstance(contigs, str):
            contigs = [contigs]
        else:
            contigs = list(contigs)
        if PAIRED and ext != "bam":
            return expand(
                os.path.join(outdir, "{{sample}}", "{{fasta}}", "{contig}{{p}}.{ext}"),
                contig=contigs,
                ext=ext,
            )
        else:
            return expand(
                os.path.join(outdir, "{{sample}}", "{{fasta}}", "{contig}.{ext}"),
                contig=contigs,
                ext=ext,
            )
    if level == "fasta":
        table = checkpoints.calculate_coverage.get(**wildcards).output[0]
        df = pd.read_csv(table, sep="\t", index_col=["samplename", "fasta"])
        fastas = list(set(df.loc[wildcards.sample].index))
        if PAIRED and ext != "bam":
            return expand(
                os.path.join(outdir, "{{sample}}", "{fasta}{{p}}.{ext}"),
                fasta=fastas,
                ext=ext,
            )
        else:
            return expand(
                os.path.join(outdir, "{{sample}}", "{fasta}.{ext}"),
                fasta=fastas,
                ext=ext,
            )


def get_header(fa):
    seqid = SeqIO.read(fa, "fasta").id
    return seqid.replace(">", "")


def seqkit_replace(wildcards):
    if PAIRED:
        return f"{wildcards.sample}_{{nr}}/{wildcards.p}"
    else:
        return f"{wildcards.sample}_{{nr}}"
