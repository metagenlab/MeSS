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


def get_fasta_dir():
    try:
        checkpoint_dir = os.path.dirname(
            checkpoints.download_assemblies.get(**wildcards).output[0]
        )
        fasta_dir = os.path.join(checkpoint_dir, "assemblies")

    except AttributeError:
        fasta_dir = os.path.abspath(config["fasta_dir"])
    return fasta_dir


def get_value(table, wildcards, value):
    if config["seq_tech"] == "illumina":
        df = pd.read_csv(table, sep="\t", index_col=["samplename", "fasta"])
        val = df.loc[wildcards.sample].loc[wildcards.fasta][value]
    else:
        df = pd.read_csv(
            table,
            sep="\t",
            dtype={"chunk": object},
            index_col=["samplename", "fasta", "chunk"],
        )
        val = df.loc[wildcards.sample].loc[wildcards.fasta].loc[wildcards.chunk][value]
    return val


def get_assembly_summary_path(wildcards):
    try:
        table = os.path.abspath(
            checkpoints.download_assemblies.get(**wildcards).output[0]
        )
    except AttributeError:
        table = os.path.abspath(config["asm_summary"])
    return table


# fasta_dir = get_fasta_dir()
# fasta_list = list_files(fasta_dir, "fasta.gz,fa.gz,fna.gz")
# fasta_names = [os.path.basename(fa).split(".fna.gz")[0] for fa in fasta_list]

# if config["paired"]:
#     pairs = [1, 2]
# else:
#     pairs = [1]


# sam_out = []
# if config["bam"] and config["paired"]:
#     sam_out = temp(f"{outdir}/fastq/{{sample}}/{{fasta}}.sam")
# elif config["bam"] and config["paired"] == False:
#     sam_out = temp(f"{outdir}/fastq/{{sample}}/{{fasta}}1.sam")
# if config["bam"] == False:
#     sam_out = temp(f"{outdir}/fastq/{{sample}}/{{fasta}}.txt")


def list_concat(wildcards, ext):
    table = checkpoints.calculate_coverage.get(**wildcards).output[0]
    df = pd.read_csv(table, sep="\t", index_col=["samplename", "fasta"])
    fastas = list(set(df.loc[wildcards.sample].index))
    if ext == "fastq":
        if config["seq_tech"] == "illumina":
            return expand(
                "{outdir}/fastq/{{sample}}/{fasta}{{p}}.fq.gz",
                outdir=outdir,
                fasta=fastas,
            )
        else:
            return expand(
                "{outdir}/fastq/{{sample}}/{fasta}.fq.gz",
                outdir=outdir,
                fasta=fastas,
            )
    elif ext == "bam":
        return expand(
            "{outdir}/bam/{{sample}}/{fasta}.bam",
            outdir=outdir,
            fasta=fastas,
        )


## PBSIM3 functions


def get_header(fa):
    seqid = SeqIO.read(fa, "fasta").id
    return seqid.replace(">", "")


def pbsim3_expand(wildcards, subdir, ext, file_type):
    table = checkpoints.calculate_coverage.get(**wildcards).output[0]
    df = pd.read_csv(table, sep="\t", index_col=["samplename", "fasta"])
    try:
        n = int(df.loc[wildcards.sample].loc[wildcards.fasta]["contig_count"].iloc[0])
    except AttributeError:
        n = int(df.loc[wildcards.sample].loc[wildcards.fasta]["contig_count"])
    contigs = [f"{x+1:04}" for x in range(n)]
    if file_type == "chunks":
        return expand(
            "{outdir}/{d}/{{sample}}/{{fasta}}-{chunk}_{{contig}}.{ext}",
            outdir=outdir,
            chunk=CHUNKS,
            ext=ext,
            d=subdir,
        )
    elif file_type == "contigs":
        return expand(
            "{outdir}/{d}/{{sample}}/{{fasta}}_{contig}.{ext}",
            outdir=outdir,
            contig=contigs,
            ext=ext,
            d=subdir,
        )
    elif file_type == "both":
        return expand(
            "{outdir}/{d}/{{sample}}/{{fasta}}-{chunk}_{contig}.{ext}",
            outdir=outdir,
            chunk=CHUNKS,
            contig=contigs,
            ext=ext,
            d=subdir,
        )
