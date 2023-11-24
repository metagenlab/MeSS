
rule sim_illumina:
    input:
        list_reads(),


rule art_illumina:
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
        -na {params.args} -o {params.prefix} &> {log}
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
