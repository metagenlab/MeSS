hifi = False
if config["seq_tech"] == "pacbio":
    model = "QSHMM-RSII"
    if config["passes"] <= 1:
        hifi = False
    elif config["passes"] >=2:
        hifi = True
    
else:
    if config["profile"] == "q20":
        model = "QSHMM-ONT-HQ"
    else:
        model = "QSHMM-ONT"

pbsim_args = f"--length-mean {config['read_len']} --length-sd {config['read_sd']}"
if hifi:
    pbsim_args += f" --pass-num {config['passes']}"

ruleorder: pbsim3 > convert_maf_to_sam

rule pbsim3:
    input:
        fa=f"{outdir}/fasta/{{fasta}}.renamed",
        df=f"{outdir}/cov.tsv",
    output:
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.sam") 
        if hifi 
        else temp(f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.fastq"),
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.maf"),
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.ref"),
    params:
        model=model,
        args=pbsim_args,
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
        pbsim --strategy wgs --method qshmm \\
        --prefix {params.prefix} --seed {params.seed} \\
        --qshmm ${{CONDA_PREFIX}}/data/{params.model}.model \\
        {params.args} --depth {params.cov} --genome {input.fa} &> {log}
        """
if hifi:
    rule ccs_sam_to_bam:
        input:
            f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.sam"
        output:
            temp(f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.bam")
        log:
            f"{outdir}/logs/ccs/{{sample}}/{{fasta}}_{{contig}}.log",

        shell:
            """
            bioconvert {input} {output} 2>> {log}
            """

    rule ccs_bam_to_fastq:
        input:
            f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.bam"
        output:
            temp(f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.fastq")
        params:
            f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.fastq.gz"
        log:
            f"{outdir}/logs/ccs/{{sample}}/{{fasta}}_{{contig}}.log",

        threads: 3
        shell:
            """
            ccs -j {threads} {input} {params}
            gunzip {params}
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
