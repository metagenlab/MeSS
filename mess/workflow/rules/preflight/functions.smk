import os
import pandas as pd
import glob
import random
from itertools import chain
from itertools import product
from Bio import SeqIO


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
    if SKIP_FA_PROC and SEQ_TECH != "illumina":
        return os.path.join(basedir, "{fasta}_{contig}.fna")


def list_fastas(wildcards, contigs=False):
    table = checkpoints.calculate_coverage.get(**wildcards).output[0]
    df = pd.read_csv(table, sep="\t")
    paths = list(set(df["path"]))
    fastas = list(set(df["fasta"]))
    if contigs:
        return expand(os.path.join(dir.out.fasta, "{fasta}.renamed"), fasta=fastas)
    else:
        return paths


def get_value(table, wildcards, value):
    df = pd.read_csv(
        table,
        sep="\t",
        index_col=["samplename", "fasta"],
    )
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


def list_cat(wildcards, ext):
    table = checkpoints.calculate_coverage.get(**wildcards).output[0]
    df = pd.read_csv(table, sep="\t", index_col=["samplename", "fasta"])
    fastas = list(set(df.loc[wildcards.sample].index))
    if ext == "fastq":
        if SEQ_TECH == "illumina":
            return expand(
                os.path.join(dir.out.short, "{{sample}}", "{fasta}{{p}}.fq.gz"),
                fasta=fastas,
            )
        else:
            return expand(
                os.path.join(dir.out.long, "{{sample}}", "{fasta}.fq.gz"),
                fasta=fastas,
            )
    elif ext == "bam":
        return expand(
            os.path.join(dir.out.bam, "{{sample}}", "{fasta}.bam"),
            fasta=fastas,
        )


## PBSIM3 functions
def pbsim3_expand(wildcards, outdir, ext):
    table = checkpoints.split_contigs.get(**wildcards).output[0]
    df = pd.read_csv(table, sep="\t", index_col="fasta")
    contigs = df.loc[wildcards.fasta]["contig"]
    if isinstance(contigs, str):
        contigs = [contigs]
    else:
        contigs = list(contigs)
    return expand(
        os.path.join(outdir, "{{sample}}", "{{fasta}}", "{contig}.{ext}"),
        outdir=outdir,
        contig=contigs,
        ext=ext,
    )


def get_header(fa):
    seqid = SeqIO.read(fa, "fasta").id
    return seqid.replace(">", "")
