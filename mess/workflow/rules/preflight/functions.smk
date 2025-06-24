import os
import datetime
import pandas as pd
import numpy as np
import glob
from itertools import chain
from itertools import product
from Bio import SeqIO
import random
import re

# To get rid of the setlocale warning
os.environ["LC_ALL"] = "C.UTF-8"


# Concatenate Snakemake's own log file with the master log file
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    shell("cat " + current_log + " >> " + LOG)


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

    if BAM or TAX:
        bams = collect(
            os.path.join(dir.out.bam, "{sample}.{bam}"),
            sample=SAMPLES,
            bam=["bam", "bam.bai"],
        )
        reads = reads + bams
    if TAX:
        tax = collect(
            os.path.join(dir.out.tax, "{sample}_{abundance}_CAMI.txt"),
            sample=SAMPLES,
            abundance=["seq", "tax"],
        )
        reads = reads + bams + tax

    if ERRFREE:
        bams_ef = expand(
            os.path.join(dir.out.ef, "{sample}.{bam}"),
            sample=SAMPLES,
            bam=["bam", "bam.bai"],
        )
        reads = reads + bams_ef

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


def strip_fasta_ext(filename):
    return re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", filename)


def get_fasta_table(wildcards):
    if "fasta_table" not in fasta_cache:
        if FASTA_DIR:
            fastas = []
            for ext in ("fa", "fasta", "fna"):
                fastas.extend(glob.glob(os.path.join(FASTA_DIR, f"*.{ext}*")))
            fa_df = pd.DataFrame([{"path": os.path.abspath(fa)} for fa in fastas])

        elif FASTA_PATH:
            fa_df = pd.read_csv(INPUT, sep="\t")[["path"]]
        else:
            fa_df = pd.read_csv(
                checkpoints.download_assemblies.get(**wildcards).output[0], sep="\t"
            )
        fa_df["fasta"] = [
            strip_fasta_ext(os.path.basename(path)) for path in fa_df["path"]
        ]
        if FASTA_DIR and "fasta" in tsv_df.columns:
            tsv_df["fasta"] = [strip_fasta_ext(fa) for fa in tsv_df["fasta"]]
            fa_df = fa_df[fa_df["fasta"].isin(tsv_df["fasta"].drop_duplicates())]

        fa_df.set_index("fasta", inplace=True)
        fasta_cache["fasta_table"] = fa_df
    fa_df = fasta_cache["fasta_table"]
    return fa_df


def fasta_input(wildcards):
    df = get_fasta_table(wildcards)
    try:
        return df.loc[wildcards.fasta]["path"].drop_duplicates()
    except AttributeError:
        return df.loc[wildcards.fasta]["path"]


def list_fastas(wildcards):
    df = get_fasta_table(wildcards)
    if PRIMERSEARCH:
        return expand(
            os.path.join(dir.out.processing, "{fasta}.amplicons.fasta"),
            fasta=list(set(df.index)),
        )
    else:
        return expand(
            os.path.join(dir.out.processing, "{fasta}.fasta"), fasta=list(set(df.index))
        )


def get_tax_date(path):
    timestamp = os.path.getctime(path)
    creation_time = datetime.datetime.fromtimestamp(timestamp)
    formatted_date = creation_time.strftime("%d-%m-%Y")
    if CUSTOM_TAX:
        return f"custom-taxonomy_{formatted_date}"
    else:
        return f"ncbi-taxonomy_{formatted_date}"


table_cache = {}


def get_cov_df(wildcards):
    return checkpoints.split_contigs.get(**wildcards).output[0]


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
    if PRIMERSEARCH or FASTA_DIR or FASTA_PATH:
        return os.path.join(dir.out.processing, "seqkit_stats.tsv")
    else:
        try:
            return checkpoints.download_assemblies.get(**wildcards).output[0]
        except AttributeError:
            return os.path.join(dir.out.processing, "seqkit_stats.tsv")


tsv_cache = {}
if "input_tsv" not in tsv_cache:
    if os.path.isfile(config.args.input):
        files = [config.args.input]
    else:
        files = glob.glob(os.path.join(config.args.input, "*.tsv"))
    tsv_df = pd.concat([pd.read_csv(file, sep="\t") for file in files])
    tsv_cache["input_tsv"] = tsv_df
else:
    tsv_df = tsv_cache["input_tsv"]

if config.args.custom_tax:
    if "custom_tax_df" not in tsv_cache:
        custom_tax_df = pd.read_csv(config.args.custom_tax, sep="\t")
        tsv_cache["custom_tax_df"] = custom_tax_df
    else:
        custom_tax_df = tsv_cache["custom_tax_df"]


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
