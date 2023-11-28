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
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.fastq"),
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.maf"),
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.ref"),
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
        f"{outdir}/logs/pbsim3/{{sample}}/{{fasta}}_{{contig}}.log",
    shell:
        """
        pbsim --strategy wgs --method errhmm \\
        --prefix {params.prefix} --seed {params.seed} \\
        --errhmm ${{CONDA_PREFIX}}/data/{params.model}.model \\
        --length-mean {params.mean_len} --length-sd {params.len_sd} \\
        --depth {params.cov} --genome {input.fa} &> {log}
        """


rule rename_maf_ref:
    input:
        maf=f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.maf",
        fa=f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.ref",
    output:
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.header"),
    params:
        seqname=lambda wildcards, input: get_header(input.fa),
    shell:
        """
        sed 's/ref/{params.seqname}/g' {input.maf} > {output}
        """


rule concat_maf:
    input:
        lambda wildcards: list_mafs(wildcards, "header"),
    output:
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}.maf"),
    shell:
        "cat {input} > {output}"


rule cat_contig_fastq:
    input:
        lambda wildcards: list_mafs(wildcards, "fastq"),
    output:
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}.fq"),
    shell:
        "cat {input} > {output}"
