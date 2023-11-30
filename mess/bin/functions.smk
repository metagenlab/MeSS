import os
import pandas as pd
import glob
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


def get_value(table, sample, fasta, value):
    df = pd.read_csv(table, sep="\t", index_col=["samplename", "fasta"])
    return df.loc[sample].loc[fasta][value]


def list_reads():
    if config["seq_tech"] == "illumina":
        reads = expand(
            "{outdir}/fastq/{sample}_R{p}.fq.gz",
            outdir=outdir,
            sample=renamed_samples,
            p=pairs,
        )
    else:
        reads = expand(
            "{outdir}/fastq/{sample}.fq.gz",
            outdir=outdir,
            sample=renamed_samples,
        )

    if config["bam"]:
        bams = expand(
            "{outdir}/bam/{sample}.{ext}",
            outdir=outdir,
            sample=renamed_samples,
            ext=["bam", "bam.bai"],
        )
        reads.append(bams)
    return reads


outdir = config["outdir"]

replicates = list(range(1, config["replicates"] + 1))

samples = []
if os.path.isfile(config["input_path"]):
    files = config["input_path"]
else:
    files = list_files(config["input_path"], "tsv")

for file in files:
    df = pd.read_csv(file, sep="\t")
    try:
        samples.append(list(set(df["sample"])))
    except KeyError:
        samples.append([os.path.basename(file).split(".")[0]])

samples = list(chain.from_iterable(samples))
samples_rep = list(product(samples, replicates))
renamed_samples = [
    "-".join([os.path.basename(sample[0]).split(".")[0], str(sample[1])])
    for sample in samples_rep
]


def get_assembly_summary_path(wildcards):
    try:
        table = os.path.abspath(
            checkpoints.download_assemblies.get(**wildcards).output[0]
        )
    except AttributeError:
        table = os.path.abspath(config["asm_summary"])
    return table


fasta_dir = get_fasta_dir()
fasta_list = list_files(fasta_dir, "fasta.gz,fa.gz,fna.gz")
fasta_names = [os.path.basename(fa).split(".fna.gz")[0] for fa in fasta_list]

if config["paired"]:
    pairs = [1, 2]
else:
    pairs = [1]


sam_out = []
if config["bam"] and config["paired"]:
    sam_out = temp(f"{outdir}/fastq/{{sample}}/{{fasta}}.sam")
elif config["bam"] and config["paired"] == False:
    sam_out = temp(f"{outdir}/fastq/{{sample}}/{{fasta}}1.sam")
if config["bam"] == False:
    sam_out = temp(f"{outdir}/fastq/{{sample}}/{{fasta}}.txt")


def list_concat(wildcards, ext):
    table = checkpoints.calculate_coverage.get(**wildcards).output[0]
    df = pd.read_csv(table, sep="\t", index_col=["samplename", "fasta"])
    fastas = list(df.loc[wildcards.sample].index)
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


# Distribute amount of passes on

if config["passes"] <= 1:
    passes_per_repeat = {1: 1}
else:
    min_passes = 2
    total_passes = config["passes"]
    repeats = list(range(1, total_passes // min_passes + 1))
    passes_per_repeat = {}
    count = total_passes
    for i in repeats:
        count = count - min_passes
        passes_per_repeat[i] = min_passes
        if (count <= min_passes) and (count > 0):
            passes_per_repeat[i] = count
        else:
            passes_per_repeat[i] = min_passes


def pbsim3_expand(wildcards, ext, file_type):
    table = checkpoints.calculate_coverage.get(**wildcards).output[0]
    df = pd.read_csv(table, sep="\t", index_col=["samplename", "fasta"])
    n = int(df.loc[wildcards.sample].loc[wildcards.fasta]["contig_count"])
    contigs = [f"{x+1:04}" for x in range(n)]
    if "fq" in ext:
        d = "fastq"
    elif "bam" in ext:
        d = "bam"
    if file_type == "passes":
        return expand(
            "{outdir}/{d}/{{sample}}/{{fasta}}-{passnum}_{{contig}}.{ext}",
            outdir=outdir,
            passnum=list(passes_per_repeat.keys()),
            ext=ext,
            d=d,
        )
    elif file_type == "contigs":
        return expand(
            "{outdir}/{d}/{{sample}}/{{fasta}}_{contig}.{ext}",
            outdir=outdir,
            contig=contigs,
            ext=ext,
            d=d,
        )
    elif file_type == "both":
        return expand(
            "{outdir}/{d}/{{sample}}/{{fasta}}-{passnum}_{contig}.{ext}",
            outdir=outdir,
            passnum=list(passes_per_repeat.keys()),
            contig=contigs,
            ext=ext,
            d=d,
        )
