
include: "./functions.smk"


rule pbsim2:
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
