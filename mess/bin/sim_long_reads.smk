if config["seq_tech"] == "pacbio":
    if config["profile"] == "sequel":
        model = "ERRHMM-SEQUEL"
    else:
        model = "ERRHMM-RSII"
else:
    if config["profile"] == "q20":
        model = "ERRHMM-ONT-HQ"
    else:
        model = "ERRHMM-ONT"


rule pbsim3:
    input:
        fa=f"{outdir}/fasta/{{fasta}}.renamed",
        df=f"{outdir}/cov.tsv",
    output:
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}_0001.fastq"),
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}_0001.maf"),
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}_0001.ref"),
    params:
        model=model,
        mean_len=config["read_len"],
        len_sd=config["read_sd"],
        cov=lambda wildcards, input: get_value(
            input.df, wildcards.sample, wildcards.fasta, "cov_sim"
        ),
        seed=lambda wildcards, input: get_value(
            input.df, wildcards.sample, wildcards.fasta, "seed"
        ),
        prefix=f"{outdir}/fastq/{{sample}}/{{fasta}}",
    log:
        f"{outdir}/logs/pbsim3/{{sample}}/{{fasta}}.log",
    shell:
        """
        pbsim --strategy wgs --method errhmm --seed {params.seed} \\
        --prefix {params.prefix} --id-prefix {wildcards.fasta} \\
        --errhmm ${{CONDA_PREFIX}}/data/{params.model}.model \\
        --length-mean {params.mean_len} --length-sd {params.len_sd} \\
        --depth {params.cov} --genome {input.fa} &> {log}
        """
