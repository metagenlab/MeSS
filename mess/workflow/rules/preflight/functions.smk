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
        files = indir
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
    products = list(product(samples, replicates))
    return [prod[0] + "-" + str(prod[1]) for prod in products]


def get_fasta(wildcards):
    try:
        checkpoint_dir = os.path.dirname(
            checkpoints.download_assemblies.get(**wildcards).output[0]
        )
        fasta_dir = os.path.asbpath(os.path.join(checkpoint_dir, "assemblies"))
        ext = "fna"

    except AttributeError:
        fasta_dir = os.path.abspath(config.args.fasta)
        fastas = list_files(fasta_dir, "gz")
        ext = list(set([os.path.basename(fa).split(".")[-2] for fa in fastas]))
        if len(ext) > 1:
            logging.error("Fastas do not have the same extensions")
        else:
            ext = ext[0]
    return os.path.join(fasta_dir, f"{{fasta}}.{ext}.gz")


def get_value(table, wildcards, value):
    if SEQ_TECH == "illumina":
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
        table = os.path.abspath(ASM_SUMMARY)
    return table


def list_concat(wildcards, ext):
    table = checkpoints.calculate_coverage.get(**wildcards).output[0]
    df = pd.read_csv(table, sep="\t", index_col=["samplename", "fasta"])
    fastas = list(set(df.loc[wildcards.sample].index))
    if ext == "fastq":
        if SEQ_TECH == "illumina":
            return expand(
                os.path.join(dir.out.fastq, "{{sample}}", "{fasta}{{p}}.fq.gz"),
                fasta=fastas,
            )
        else:
            return expand(
                os.path.join(dir.out.fastq, "{{sample}}", "{fasta}.fq.gz"),
                fasta=fastas,
            )
    elif ext == "bam":
        return expand(
            os.path.join(dir.out.fastq, "{{sample}}", "{fasta}.bam"),
            fasta=fastas,
        )


## PBSIM3 functions


def get_header(fa):
    seqid = SeqIO.read(fa, "fasta").id
    return seqid.replace(">", "")


def pbsim3_expand(wildcards, subdir, ext):
    table = checkpoints.calculate_coverage.get(**wildcards).output[0]
    df = pd.read_csv(table, sep="\t", index_col=["samplename", "fasta"])
    try:
        n = int(df.loc[wildcards.sample].loc[wildcards.fasta]["contig_count"].iloc[0])
    except AttributeError:
        n = int(df.loc[wildcards.sample].loc[wildcards.fasta]["contig_count"])
    contigs = [f"{x+1:04}" for x in range(n)]
    return expand(
        os.path.join(
            dir.out.base, ext, "{{sample}}", "{{fasta}}-{chunk}_{contig}.{ext}"
        ),
        chunk=CHUNKS,
        contig=contigs,
        ext=ext,
        subdir=subdir,
    )
