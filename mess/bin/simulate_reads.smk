import pandas as pd
import os
import random
import glob

random.seed(config["seed"])
replicates = list(range(1, config["replicates"] + 1))
community_name = config["community_name"]
read_direction = []
if config["seq_tech"] != "illumina":
    read_direction = config["seq_tech"]
if config["seq_tech"] == "illumina" and config["read_status"] == "single":
    read_direction = "single"
if config["seq_tech"] == "illumina" and config["read_status"] == "paired":
    read_direction = ["R1", "R2"]
pipeline_path = workflow.basedir
seed_list = random.sample(range(1, 1000000), len(replicates))
seed_dic = dict(zip(replicates, seed_list))


def get_input_value(path):
    """
    Function that parses input table and raises exceptions if there are NA in the table
    Returns the input value set by the user, useful for read count calculation
    """
    tb = pd.read_csv(path, sep="\t")
    acceptedvals = ["Reads", "Coverage", "ReadPercent", "RelativeProp"]
    lastcol = tb.columns[len(tb.columns) - 1]
    nans = tb[lastcol].isna().sum()
    if nans > 0:
        raise Exception(f"{nans} values in {lastcol} are NA, make sure to input values")
    try:
        lastcol = config[
            "dist"
        ]  # If the user has set a specific distribution for read counts
        return lastcol
    except KeyError:
        if lastcol in acceptedvals:
            return lastcol
        else:
            raise Exception(
                f"{lastcol} is not recognized as an input to simulate reads"
            )


def get_read_counts_seed(rep, rep2seed, general_seed):
    """
    Function that returns a seed for calculating read counts
    By default each replicate has a unique seed
    If the param same_dist_rep is True, return a single seed for all replicates
    """
    try:
        same_dist_rep = config["same_dist_rep"]
        if same_dist_rep:
            return general_seed
    except KeyError:
        return rep2seed[int(rep)]


input_val = get_input_value(config["input_table_path"])


rule get_entries:
    input:
        config["input_table_path"],
    output:
        temp(f"{community_name}/entries.tsv"),
    run:
        pd.read_csv(input[0], sep="\t").iloc[:, 0].to_csv(
            output[0], sep="\t", index=None
        )


checkpoint download_assemblies:
    input:
        f"{community_name}/entries.tsv",
    output:
        f"{community_name}/assembly_summary.tsv",
    params:
        ncbi_key=config["NCBI_key"],
        ncbi_email=config["NCBI_email"],
        prefix=f"{community_name}/download",
    shell:
        """
        assembly_finder -i {input} -o {params.prefix} -nk {params.ncbi_key} -ne {params.ncbi_email} nolock
        """


rule create_read_counts_table:
    """
    Rule for generating a table per sample with read counts for each genome 
    """
    input:
        input_tb=config["input_table_path"],
        assemblies_tb=f"{community_name}/assembly_summary.tsv",
    output:
        rc_table="readcounts-{community}-{rep}.tsv",
    params:
        value=input_val,
        seed=lambda wildcards: get_read_counts_seed(
            wildcards.rep, seed_dic, config["seed"]
        ),
    log:
        "logs/read_counts_table/{community}-{rep}.txt",
    script:
        "read_counts_table.py"


rule decompress_assemblies:
    input:
        f"{community_name}/assemblies/{{assemblyname}}_genomic.fna.gz",
    output:
        temp(
            f"{community_name}/assemblies/{{assemblyname,[0-9a-zA-Z._-]+}}_genomic.fna"
        ),
    benchmark:
        "benchmark/decompress/{assemblyname}.txt"
    shell:
        "zcat {input[0]} > {output}"


rule merge_contigs:
    input:
        f"{community_name}/assemblies/{{assemblyname}}_genomic.fna",
    output:
        temp("{community,[-_0-9a-zA-Z]+}-{assemblyname,[0-9a-zA-Z._-]+}.fa"),
    benchmark:
        "benchmark/merge/{community}-{assemblyname}.txt"
    script:
        "merge_contigs.py"


if config["seq_tech"] == "illumina":

    def get_read_num(wildcards, table):
        if not os.path.exists(table):
            return -1
        with open(table) as f:
            tb = pd.read_csv(f, delimiter="\t", index_col="AssemblyNames")
            reads = tb.loc[wildcards.assemblyname]["Reads"]
            return reads

    if config["read_status"] == "paired":

        rule generate_illumina_paired_reads:
            input:
                fa="{community}-{assemblyname}.fa",
                tab="readcounts-{community}-{rep}.tsv",
            output:
                temp(
                    "simreads/{community,[-_0-9a-zA-Z]+}-{rep,[0-9]+}-{assemblyname,[0-9a-zA-Z._-]+}_R1.fq"
                ),
                temp(
                    "simreads/{community,[-_0-9a-zA-Z]+}-{rep,[0-9]+}-{assemblyname,[0-9a-zA-Z._-]+}_R2.fq"
                ),
                temp(
                    "simreads/{community,[-_0-9a-zA-Z]+}-{rep,[0-9]+}-{assemblyname,[0-9a-zA-Z._-]+}_R.sam"
                ),
            params:
                seq_system=config["illumina_sequencing_system"],
                read_length=config["illumina_read_len"],
                mean_frag_len=config["illumina_mean_frag_len"],
                sd=config["illumina_sd_frag_len"],
                read_num=lambda wildcards, input: get_read_num(wildcards, input.tab),
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

    if config["read_status"] == "single":

        rule generate_illumina_single_reads:
            input:
                fa="assembly_gz/{community}-{assemblyname}.fa",
                tab="readcounts-{community}-{rep}.tsv",
            output:
                temp(
                    "simreads/{community,[-_0-9a-zA-Z]+}-{rep,[0-9]+}-{assemblyname,[0-9a-zA-Z._-]+}_single.fq"
                ),
                temp(
                    "simreads/{community,[-_0-9a-zA-Z]+}-{rep,[0-9]+}-{assemblyname,[0-9a-zA-Z._-]+}_single.sam"
                ),
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
            temp(
                "simreads/{community,[-_0-9a-zA-Z]+}-{rep,[0-9]+}-{assemblyname,[0-9a-zA-Z._-]+}_0001.fastq"
            ),
            temp(
                "simreads/{community,[-_0-9a-zA-Z]+}-{rep,[0-9]+}-{assemblyname,[0-9a-zA-Z._-]+}_0001.maf"
            ),
            temp(
                "simreads/{community,[-_0-9a-zA-Z]+}-{rep,[0-9]+}-{assemblyname,[0-9a-zA-Z._-]+}_0001.ref"
            ),
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
        checkpoint_directory = f"{community_name}/assemblies"
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
