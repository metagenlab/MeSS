import os
import pandas as pd
import glob
from itertools import chain
from itertools import product


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
            "{outdir}/bam/{sample}.bam", outdir=outdir, sample=renamed_samples
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

art_args = ""
if config["paired"]:
    art_args += f"-p -m {config['mean_frag_len']} -s {config['sd_frag_len']}"
if config["bam"]:
    art_args += " -sam"


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
        return expand(
            "{outdir}/fastq/{{sample}}/{fasta}{{p}}.fq",
            outdir=outdir,
            fasta=fastas,
        )
    elif ext == "bam":
        return expand(
            "{outdir}/bam/{{sample}}/{fasta}.bam",
            outdir=outdir,
            fasta=fastas,
        )
