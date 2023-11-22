include: "prepare_tables.smk"


import os
import json
import glob


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


def get_json_value(json_file, fasta, value):
    d = json.load(open(json_file))
    return d[value][fasta]


fasta_dir = get_fasta_dir()
fasta_list = list_files(fasta_dir, "fasta.gz,fa.gz,fna.gz")
fasta_names = [os.path.basename(fa).split(".fna.gz")[0] for fa in fasta_list]

if config["paired"]:
    pairs = [1, 2]
else:
    pairs = [1]


def list_illumina_reads(bam):
    reads = expand(
        "{outdir}/fastq/{sample}_R{p}.fq.gz",
        outdir=outdir,
        sample=renamed_samples,
        p=pairs,
    )
    if bam:
        bams = expand(
            "{outdir}/bam/{sample}.bam", outdir=outdir, sample=renamed_samples
        )
        reads.append(bams)
    return reads


rule sim_illumina:
    input:
        list_illumina_reads(config["bam"]),


rule decompress_assemblies:
    input:
        f"{fasta_dir}/{{fasta}}.fna.gz",
    output:
        temp(f"{outdir}/fasta/{{fasta}}.fa"),
    benchmark:
        f"benchmark/decompress/{outdir}/{{fasta}}.txt"
    shell:
        "zcat {input} > {output}"


rule merge_contigs:
    input:
        f"{outdir}/fasta/{{fasta}}.fa",
    output:
        temp(f"{outdir}/fasta/{{fasta}}.merged"),
    benchmark:
        "benchmark/merge/{fasta}.txt"
    script:
        "merge_contigs.py"


art_args = {"bam": "-sam", "paired": "-p"}


def get_art_args(args):
    return " ".join([art_args[arg] for arg in args if config[arg]])


rule sim_art_illumina_reads:
    input:
        fasta=f"{outdir}/fasta/{{fasta}}.merged",
        json=f"{outdir}/{{sample}}-cov.json",
    output:
        fastqs=[
            temp(f"{outdir}/fastq/{{sample}}/{{fasta}}1.fq.gz"),
            temp(f"{outdir}/fastq/{{sample}}/{{fasta}}2.fq.gz"),
        ]
        if config["paired"]
        else temp(f"{outdir}/fastq/{{sample}}/{{fasta}}1.fq.gz"),
        sam=temp(f"{outdir}/fastq/{{sample}}/{{fasta}}.sam")
        if config["paired"]
        else temp(f"{outdir}/fastq/{{sample}}/{{fasta}}1.sam"),
    params:
,
        args=get_art_args(["bam", "paired"]),
        seq_system=config["seq_system"],
        read_len=config["read_len"],
        cov=lambda wildcards, input: get_json_value(
            input.json, wildcards.fasta, "cov_sim"
        ),
        frag_len=config["mean_frag_len"],
        sd=config["sd_frag_len"],
        seed=lambda wildcards, input: get_json_value(
            input.json, wildcards.fasta, "seed"
        ),
        prefix=f"{outdir}/fastq/{{sample}}/{{fasta}}"
        if config["paired"]
        else f"{outdir}/fastq/{{sample}}/{{fasta}}1",
    log:
        "logs/art_illumina/{sample}/{fasta}.log",
    benchmark:
        "benchmark/art_illumina/{sample}/{fasta}.txt"
    shell:
        """
        art_illumina -ss {params.seq_system} \\
        -i {input.fasta} -rs {params.seed} \\
        -m {params.frag_len} -s {params.sd} \\
        -l {params.read_len} -f {params.cov} \\
        -d {wildcards.fasta} -sam -na \\
        -o {params.prefix} {params.args} &> {log}
        gzip {params.prefix}*.fq
        """


rule convert_sam_to_bam:
    input:
        f"{outdir}/fastq/{{sample}}/{{fasta}}.sam"
        if config["paired"]
        else f"{outdir}/fastq/{{sample}}/{{fasta}}1.sam",
    output:
        f"{outdir}/bam/{{sample}}/{{fasta}}.bam",
    log:
        "logs/bam/{sample}/{fasta}.log",
    shell:
        """
        bioconvert {input} {output} 2> {log}
        """


def list_concat(file, sample, ext):
    d = json.load(open(file))
    fastas = d[sample]
    if ext == "fastq":
        return expand(
            "{outdir}/fastq/{{sample}}/{fasta}{{p}}.fq.gz",
            outdir=outdir,
            fasta=fastas,
        )
    elif ext == "bam":
        return expand(
            "{outdir}/bam/{{sample}}/{fasta}.bam",
            outdir=outdir,
            fasta=fastas,
        )


rule concat_bam:
    input:
        lambda wildcards: list_concat(
            f"{outdir}/sample2fa.json", wildcards.sample, "bam"
        ),
    output:
        f"{outdir}/bam/{{sample}}.bam",
    shell:
        """
        samtools merge -o {output} {input}
        """


rule concat_fastq:
    input:
        lambda wildcards: list_concat(
            f"{outdir}/sample2fa.json", wildcards.sample, "fastq"
        ),
    output:
        f"{outdir}/fastq/{{sample}}_R{{p}}.fq.gz",
    shell:
        """
        cat {input} {output}
        """
