art_args = ""
if config["paired"]:
    art_args += f"-p -m {config['frag_len']} -s {config['frag_sd']}"
if config["bam"]:
    art_args += " -sam"


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
        seq_system=config["profile"],
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
