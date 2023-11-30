hifi = False
if config["passes"] > 1:
    hifi = True
if config["seq_tech"] == "pacbio":
    model = "QSHMM-RSII"

else:
    if config["profile"] == "R10.4":
        model = "QSHMM-ONT-HQ"


wildcard_constraints:
    contig="[0-9]+",
    sample="[^/]+",
    fasta="[^/-]+",
    passnum="[0-9]+",


rule copy_fasta_per_passes:
    input:
        f"{outdir}/fasta/{{fasta}}.renamed",
    output:
        temp(f"{outdir}/fasta/{{fasta}}-{{passnum}}.fa"),
    shell:
        "cp {input} {output}"


rule pbsim3:
    input:
        fa=f"{outdir}/fasta/{{fasta}}-{{passnum}}.fa",
        df=f"{outdir}/cov.tsv",
    output:
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}-{{passnum}}_{{contig}}.sam")
        if hifi
        else temp(f"{outdir}/fastq/{{sample}}/{{fasta}}-{{passnum}}_{{contig}}.fastq"),
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}-{{passnum}}_{{contig}}.maf"),
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}-{{passnum}}_{{contig}}.ref"),
    params:
        model=model,
        mean_len=config["read_len"],
        len_sd=config["read_sd"],
        pass_num=config["passes"],
        cov=lambda wildcards, input: get_value(
            input.df, wildcards.sample, wildcards.fasta, "cov_sim"
        ),
        seed=lambda wildcards, input: get_value(
            input.df, wildcards.sample, wildcards.fasta, "seed"
        ),
        prefix=f"{outdir}/fastq/{{sample}}/{{fasta}}-{{passnum}}",
    log:
        f"{outdir}/logs/pbsim3/{{sample}}/{{fasta}}-{{passnum}}_{{contig}}.log",
    shell:
        """
        pbsim --strategy wgs --method qshmm \\
        --prefix {params.prefix} --seed {params.seed} \\
        --qshmm ${{CONDA_PREFIX}}/data/{params.model}.model \\
        --length-mean {params.mean_len} --length-sd {params.len_sd} \\
        --pass-num {params.pass_num} \\
        --depth {params.cov} --genome {input.fa} &> {log}
        """


if hifi:

    rule cat_ccs_sam:
        input:
            lambda wildcards: pbsim3_expand(wildcards, "sam", "both"),
        output:
            temp(f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.sam"),
        shell:
            "cat {input} {output}"

    rule ccs_sam_to_bam:
        input:
            f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.sam",
        output:
            temp(f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.bam"),
        log:
            f"{outdir}/logs/bam/{{sample}}/{{fasta}}_{{contig}}.log",
        threads: 6
        shell:
            """
            samtools view -@ {threads} -bS {input} \
            | samtools sort -@ {threads} > {output} 2>> {log}
            """

    rule ccs_bam_to_fastq:
        input:
            f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.bam",
        output:
            temp(f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.fq.gz"),
        log:
            f"{outdir}/logs/ccs/{{sample}}/{{fasta}}_{{contig}}.log",
        threads: 6
        shell:
            """
            ccs -j {threads} {input} {output} &> {log}
            """


rule add_reference_name:
    input:
        maf=f"{outdir}/fastq/{{sample}}/{{fasta}}-{{passnum}}_{{contig}}.maf",
        fa=f"{outdir}/fastq/{{sample}}/{{fasta}}-{{passnum}}_{{contig}}.ref",
    output:
        temp(f"{outdir}/bam/{{sample}}/{{fasta}}-{{passnum}}_{{contig}}.maf"),
    params:
        seqname=lambda wildcards, input: get_header(input.fa),
    shell:
        """
        sed 's/ref/{params.seqname}/g' {input.maf} > {output}
        """


rule convert_maf_to_bam:
    input:
        f"{outdir}/bam/{{sample}}/{{fasta}}-{{passnum}}_{{contig}}.maf",
    output:
        sam=temp(f"{outdir}/bam/{{sample}}/{{fasta}}-{{passnum}}_{{contig}}.sam"),
        bam=temp(f"{outdir}/bam/{{sample}}/{{fasta}}-{{passnum}}_{{contig}}.bam"),
    log:
        f"{outdir}/logs/bam/{{sample}}/{{fasta}}-{{passnum}}_{{contig}}.log",
    threads: 3
    shell:
        """
        bioconvert {input} {output.sam} 2>> {log}
        samtools view -@ {threads} -bS {output.sam} \\
        | samtools sort -@ threads > {output.bam} 2>> {log}
        """


rule cat_pbsim_bam:
    input:
        lambda wildcards: pbsim3_expand(wildcards, "bam", "both"),
    output:
        temp(f"{outdir}/bam/{{sample}}/{{fasta}}.bam"),
    log:
        f"{outdir}/logs/bam/{{sample}}/{{fasta}}.log",
    threads: 3
    shell:
        "samtools merge -@ {threads} -o {output} {input} 2>> {log}"


if not hifi:

    rule compress_single_pass_fastq:
        input:
            f"{outdir}/fastq/{{sample}}/{{fasta}}-1_{{contig}}.fastq",
        output:
            temp(f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.fq.gz"),
        params:
            f"{outdir}/fastq/{{sample}}/{{fasta}}_{{contig}}.fq",
        threads: 3
        shell:
            """
            mv {input} {params}
            pigz -p {threads} {params}
            """


rule cat_contigs_fastq:
    input:
        lambda wildcards: pbsim3_expand(wildcards, "fq.gz", "contigs"),
    output:
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}.fq.gz"),
    shell:
        "cat {input} > {output}"
