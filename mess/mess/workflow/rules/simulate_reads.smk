import json
import pandas as pd
import os
import random
import glob
from itertools import chain
from itertools import product


def list_files(indir, extensions):
    extensions = extensions.split(",")
    files = [glob.glob(f"{indir}/*.{e}") for e in extensions]
    return list(chain.from_iterable(files))


def get_cov_seed(rep, sample2seed, general_seed):
    """
    Function that returns a seed for coverage
    By default each replicate has a unique seed
    If the param same_dist_rep is True, return a single seed for all replicates
    """
    try:
        same_dist_rep = config["same_dist_rep"]
        if same_dist_rep:
            return general_seed
    except KeyError:
        return sample2seed[int(rep)]


outdir = config["outdir"]
random.seed(config["seed"])
replicates = list(range(1, config["replicates"] + 1))
samples = list_files(config["input_path"], "tsv")
sample_list = list(product(samples, replicates))

seed_list = random.sample(range(1, 1000000), len(sample_list) * len(replicates))
seed_dic = dict(zip(replicates, seed_list))


checkpoint download_assemblies:
    input:
        config["input_path"],
    output:
        f"{outdir}/download/assembly_summary.tsv",
    params:
        ncbi_key=config["NCBI_key"],
        ncbi_email=config["NCBI_email"],
        prefix=f"{outdir}/download",
    shell:
        """
        assembly_finder -i {input} -o {params.prefix} -nk {params.ncbi_key} -ne {params.ncbi_email} --nolock
        """


def get_assembly_summary_path(wildcards):
    try:
        table = os.path.abspath(
            checkpoints.download_assemblies.get(**wildcards).output[0]
        )
    except KeyError:
        table = os.path.abspath(config["fasta_dir"])
    return table


rule calculate_coverage:
    input:
        entry=lambda wildcards: list_files(config["input_path"], ".tsv"),
        asm=lambda wildcards: get_assembly_summary_path(wildcards),
    output:
        f"{outdir}/{{sample}}-{{rep}}-cov.json",
    params:
        dist=config["dist"],
        mu=config["mu"],
        sigma=config["sigma"],
        total_bases=config["total_bases"],
        read_len=config["read_len"],
        pairing=config["paired"],
        rep_sd=config["rep_sd"],
        seed=lambda wildcards: get_cov_seed(wildcards.rep, seed_dic, config["seed"]),
    log:
        f"logs/tables/{{sample}}-{{rep}}.tsv",
    shell:
        """
        calculate_cov.py -i {input.entry} -a {input.asm} \
        -d {params.dist} -m {params.mu} -s {params.sigma} \
        -n {params.total_bases} -l {params.read_len} -p {params.pairing} \
        -sd {params.rep_sd} --seed {params.seed}
        """


def get_fasta_dir(wildcards):
    try:
        table_dir = os.path.dirname(
            checkpoints.download_assemblies.get(**wildcards).output[0]
        )
        fasta_dir = os.path.join(table_dir, "assemblies")
    except KeyError:
        fasta_dir = os.path.abspath(config["fasta_dir"])
    return fasta_dir


fasta_dir = lambda wildcards: get_fasta_dir(wildcards)


rule decompress_assemblies:
    input:
        f"{fasta_dir}/{{fasta}}.gz",
    output:
        temp(f"{outdir}/{{fasta}}"),
    benchmark:
        f"benchmark/decompress/{outdir}/{{fasta}}.txt"
    shell:
        "zcat {input} > {output}"


def get_coverage(json, fasta):
    d = json.load(open(json))
    return d["cov_sim"][fasta]


if config["seq_tech"] == "illumina":

    rule merge_contigs:
        input:
            f"{outdir}/{{fasta}}",
        output:
            temp(f"{outdir}/{{fasta}}.merged"),
        benchmark:
            "benchmark/merge/{fasta}.txt"
        script:
            "merge_contigs.py"

    if config["paired"]:

        rule generate_illumina_paired_reads:
            input:
                fa=f"{outdir}/{{fasta}}.merged",
                json=f"{outdir}/coverage-{{rep}}.json",
            output:
                temp(f"{outdir}/fastq/{{fasta}}_R1.fq"),
                temp(f"{outdir}/fastq/{{fasta}}_R2.fq"),
                temp(f"{outdir}/fastq/{{fasta}}_R.sam"),
            params:
                seq_system=config["illumina_sequencing_system"],
                read_length=config["read_len"],
                mean_frag_len=config["illumina_mean_frag_len"],
                sd=config["illumina_sd_frag_len"],
                cov=lambda wildcards, input: get_coverage(input.json, wildcards.fasta),
                random_seed=lambda wildcards: seed_dic[int(wildcards.rep)],
            resources:
                nb_simulation=1,
            log:
                "logs/art_illumina/{community}-replicate-{rep}-{assemblyname}.log",
            benchmark:
                "benchmark/art_illumina/{community}-{rep}-{assemblyname}.txt"
            shell:
                """
                art_illumina -d {wildcards.assemblyname} -rs {params.random_seed} -ss {params.seq_system} -na -sam \
                -i {input.fa} -p -nf 1 -m {params.mean_frag_len} -s {params.sd}  -l {params.read_length} -c {params.read_num} \
                -o simreads/{wildcards.community}-{wildcards.rep}-{wildcards.assemblyname}_R &> {log}
                """

    else:

        rule generate_illumina_single_reads:
            input:
                fa="assembly_gz/{community}-{assemblyname}.fa",
                tab="readcounts-{community}-{rep}.tsv",
            output:
                temp("{outdir}/fastq/{{fasta}}_single.fq"),
                temp("{outdir}/fastq/{{fasta}}_single.sam"),
            params:
                seq_system=config["illumina_sequencing_system"],
                read_length=config["illumina_read_len"],
                read_num=lambda wildcards, input: get_read_num(wildcards, input.tab),
                random_seed=lambda wildcards: seed_dic[int(wildcards.rep)],
            resources:
                nb_simulation=1,
            log:
                "logs/art_illumina/{community}-{rep}-{assemblyname}.log",
            benchmark:
                "benchmark/art_illumina/{community}-{rep}-{assemblyname}.txt"
            shell:
                """
                art_illumina -d {wildcards.assemblyname} -rs {params.random_seed} -ss {params.seq_system} -sam \
                -na -i {input.fa} -nf 1 -l {params.read_length} -c {params.read_num} \
                -o simreads/{wildcards.community}-{wildcards.rep}-{wildcards.assemblyname}_single &> {log}
                """

    def get_sam(pairing):
        if pairing == "single":
            return "simreads/{community}-{rep}-{assemblyname}_single.sam"
        else:
            return "simreads/{community}-{rep}-{assemblyname}_R.sam"

    rule convert_sam_to_bam_files:
        input:
            get_sam(config["read_status"]),
        output:
            "bam/{community,[-_0-9a-zA-Z]+}-{rep,[0-9]+}-{assemblyname,[0-9a-zA-Z._-]+}.bam",
        log:
            "logs/bam/{community}-{rep}-{assemblyname}.log",
        shell:
            """
            bioconvert {input} {output} 2> {log}
            """

elif config["seq_tech"] == "ont" or config["seq_tech"] == "pacbio":

    def get_coverage(wildcards, table):
        if not os.path.exists(table):
            return -1
        with open(table) as f:
            tb = pd.read_csv(f, delimiter="\t", index_col="AssemblyNames")
            cov = tb.loc[wildcards.assemblyname]["Coverage"]
            return cov

    rule generate_long_reads:
        input:
            fa="assembly_gz/{community}-{assemblyname}.fa",
            tab="readcounts-{community}-{rep}.tsv",
        output:
            temp("{outdir}/fastq/{{fasta}}_0001.fastq"),
            temp("{outdir}/fastq/{{fasta}}_0001.maf"),
            temp("{outdir}/fastq/{{fasta}}_0001.ref"),
        params:
            chemistry=config["chemistry"],
            min_read_len=config["longreads_min_len"],
            max_read_len=config["longreads_max_len"],
            sd_read_len=config["longreads_sd_len"],
            mean_read_len=config["longreads_mean_len"],
            accuracy=config["longreads_mean_acc"],
            diff_ratio=config["difference_ratio"],
            coverage=lambda wildcards, input: get_coverage(wildcards, input.tab),
            seed=lambda wildcards: seed_dic[int(wildcards.rep)],
        resources:
            nb_simulation=1,
        log:
            "logs/pbsim2/{community}-replicate-{rep}-{assemblyname}.log",
        benchmark:
            "benchmark/pbsim2/{community}-{rep}-{assemblyname}.txt"
        shell:
            """
            pbsim --prefix simreads/{wildcards.community}-{wildcards.rep}-{wildcards.assemblyname} \
            --id-prefix {wildcards.assemblyname} --depth {params.coverage} --length-min {params.min_read_len} \
            --length-max {params.max_read_len} --difference-ratio {params.diff_ratio} --seed {params.seed} \
            --hmm_model ${{CONDA_PREFIX}}/data/{params.chemistry}.model --length-mean {params.mean_read_len} \
            --length-sd {params.sd_read_len} {input.fa} &> {log}
            """

    rule convert_maf_to_bam_files:
        input:
            "simreads/{community}-{rep}-{assemblyname}_0001.maf",
        output:
            sam=temp(
                "bam/{community,[-_0-9a-zA-Z]+}-{rep,[0-9]+}-{assemblyname,[0-9a-zA-Z._-]+}.sam"
            ),
            bam="bam/{community,[-_0-9a-zA-Z]+}-{rep,[0-9]+}-{assemblyname,[0-9a-zA-Z._-]+}.bam",
        log:
            "logs/bam/{community}-{rep}-{assemblyname}.log",
        shell:
            """
            bioconvert {input} {output.sam} 2>> {log}
            bioconvert {output.sam} {output.bam} 2>> {log}
            """


def concat_fastq(wildcards):
    try:
        directory = config["assemblies_dir"]
        files = glob.glob(os.path.join(directory, "*_genomic.fna*"))
        assemblynames = [file.split("/")[-1].split("_genomic.fna")[0] for file in files]
    except KeyError:
        checkpoint_directory = f"{fasta_dir}"
        assemblynames = glob_wildcards(
            os.path.join(checkpoint_directory, "{i}_genomic.fna.gz")
        ).i
        print("assemblynames", assemblynames)
    if config["seq_tech"] == "illumina":
        return expand(
            "simreads/{community}-{{rep}}-{assemblyname}_{{rd}}.fq",
            community=community_name,
            assemblyname=assemblynames,
        )
    elif config["seq_tech"] == "ont" or config["seq_tech"] == "pacbio":
        return expand(
            "simreads/{community}-{{rep}}-{assemblyname}_0001.fastq",
            community=community_name,
            assemblyname=assemblynames,
        )


rule concat_fastqs:
    input:
        concat_fastq,
    output:
        temp("simreads/{community}-{rep}_{rd}.unshuffled.fq"),
    resources:
        parallel_cat=1,
    benchmark:
        "benchmark/concat/{community}-{rep}_{rd}.txt"
    shell:
        """
        cat {input} > {output}
        """


rule shuffle_fastqs:
    input:
        "simreads/{community}-{rep}_{rd}.unshuffled.fq",
    output:
        "simreads/{community}-{rep}_{rd}.fq.gz",
    params:
        seed=lambda wildcards: seed_dic[int(wildcards.rep)],
    log:
        "logs/shuffling/{community}-{rep}_{rd}.log",
    benchmark:
        "benchmark/shuffle/{community}-{rep}_{rd}.txt"
    threads: 10
    shell:
        """
        seqkit shuffle {input} -j {threads} -s {params.seed} -o {output}  &> {log}
        """


rule get_taxonomic_profile:
    input:
        reads="simreads/{community}-{rep}_{rd}.fq.gz",
        table="readcounts-{community}-{rep}.tsv",
    output:
        tax="taxonomy-{community}-{rep}_{rd}.tsv",
        krona="krona/{community}-{rep}_{rd}.txt",
    params:
        config["seq_tech"],
    benchmark:
        "benchmark/taxonomic_profile/{community}-{rep}_{rd}.txt"
    script:
        "true-read-counts.py"


rule create_krona_chart:
    input:
        "krona/{community}-{rep}_{rd}.txt",
    output:
        "krona/{community}-{rep}_{rd}.html",
    log:
        "logs/krona/{community}-{rep}_{rd}.log",
    shell:
        "ktImportText -o {output} {input} &> {log}"
