include: "prepare_tables.smk"


import os
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


def get_value(table, sample, fasta, value):
    df = pd.read_csv(table, sep="\t", index_col=["samplename", "fasta"])
    return df.loc[sample].loc[fasta][value]


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
    shell:
        "zcat {input} > {output}"


rule merge_contigs:
    input:
        f"{outdir}/fasta/{{fasta}}.fa",
    output:
        temp(f"{outdir}/fasta/{{fasta}}.merged"),
    script:
        "merge_contigs.py"


art_args = {"bam": "-sam", "paired": "-p"}


def get_art_args(args):
    return " ".join([art_args[arg] for arg in args if config[arg]])


print(get_art_args)


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


rule sim_art_illumina_reads:
    input:
        fasta=f"{outdir}/fasta/{{fasta}}.merged",
        df=f"{outdir}/cov.tsv",
    output:
        sam=sam_out,
        fastqs=[
            temp(f"{outdir}/fastq/{{sample}}/{{fasta}}1.fq"),
            temp(f"{outdir}/fastq/{{sample}}/{{fasta}}2.fq"),
        ]
        if config["paired"]
        else temp(f"{outdir}/fastq/{{sample}}/{{fasta}}1.fq"),
    params:
        args=art_args,
        seq_system=config["seq_system"],
        read_len=config["read_len"],
        cov=lambda wildcards, input: get_value(
            input.df, wildcards.sample, wildcards.fasta, "cov_sim"
        ),
        seed=lambda wildcards, input: get_value(
            input.df, wildcards.sample, wildcards.fasta, "seed"
        ),
        prefix=f"{outdir}/fastq/{{sample}}/{{fasta}}"
        if config["paired"]
        else f"{outdir}/fastq/{{sample}}/{{fasta}}1",
    log:
        "logs/art_illumina/{sample}/{fasta}.log",
    shell:
        """
        art_illumina -ss {params.seq_system} \\
        -i {input.fasta} -rs {params.seed} \\
        -l {params.read_len} -f {params.cov} \\
        -d {wildcards.fasta} -na \\
        -o {params.prefix} {params.args} &> {log}
        touch {output.sam}  
        """


rule convert_sam_to_bam:
    input:
        f"{outdir}/fastq/{{sample}}/{{fasta}}.sam"
        if config["paired"]
        else f"{outdir}/fastq/{{sample}}/{{fasta}}1.sam",
    output:
        temp(f"{outdir}/bam/{{sample}}/{{fasta}}.bam"),
    log:
        "logs/bam/{sample}/{fasta}.log",
    shell:
        """
        bioconvert {input} {output} 2> {log}
        """


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


rule concat_bam:
    input:
        lambda wildcards: list_concat(wildcards, "bam"),
    output:
        f"{outdir}/bam/{{sample}}.bam",
    threads: 3
    shell:
        """
        samtools merge -@ {threads} -o {output} {input}
        """


rule concat_fastq:
    input:
        lambda wildcards: list_concat(wildcards, "fastq"),
    output:
        f"{outdir}/fastq/{{sample}}_R{{p}}.fq.gz",
    threads: 3
    shell:
        """
        cat {input} | pigz -p {threads} > {output}
        """
