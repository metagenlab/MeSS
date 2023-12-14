if config["seq_tech"] == "illumina":

    rule convert_sam_to_bam:
        input:
            f"{outdir}/fastq/{{sample}}/{{fasta}}1.sam"
            if (not config["paired"])
            else f"{outdir}/fastq/{{sample}}/{{fasta}}.sam",
        output:
            temp(f"{outdir}/bam/{{sample}}/{{fasta}}.bam"),
        log:
            f"{outdir}/logs/bam/{{sample}}/{{fasta}}.log",
        threads: 3
        shell:
            """
            samtools view -@ {threads} -bS {input} | \
            samtools sort -@ {threads} > {output} 2> {log}
            """

    rule compress_fastq:
        input:
            fastqs=[
                temp(f"{outdir}/fastq/{{sample}}/{{fasta}}1.fq"),
                temp(f"{outdir}/fastq/{{sample}}/{{fasta}}2.fq"),
            ]
            if config["paired"]
            else temp(f"{outdir}/fastq/{{sample}}/{{fasta}}1.fq"),
        output:
            fastqs=[
                temp(f"{outdir}/fastq/{{sample}}/{{fasta}}1.fq.gz"),
                temp(f"{outdir}/fastq/{{sample}}/{{fasta}}2.fq.gz"),
            ]
            if config["paired"]
            else temp(f"{outdir}/fastq/{{sample}}/{{fasta}}1.fq"),
        threads: 3
        shell:
            """
            pigz -p {threads} *.fq
            """


rule concat_bam:
    input:
        lambda wildcards: list_concat(wildcards, "bam"),
    output:
        f"{outdir}/bam/{{sample}}.bam",
    threads: 3
    log:
        f"{outdir}/logs/bam/{{sample}}.log",
    shell:
        """
        samtools merge -@ {threads} -o {output} {input} 2>> {log}
        """


rule index_bam:
    input:
        f"{outdir}/bam/{{sample}}.bam",
    output:
        f"{outdir}/bam/{{sample}}.bam.bai",
    threads: 3
    shell:
        """
        samtools index -@ {threads} {input}
        """


rule cat_fastqs:
    input:
        lambda wildcards: list_concat(wildcards, "fastq"),
    output:
        f"{outdir}/unshuffled/{{sample}}_R{{p}}.fq.gz"
        if config["seq_tech"] == "illumina"
        else f"{outdir}/unshuffled/{{sample}}.fq.gz",
    shell:
        """
        cat {input} > {output}
        """


rule shuffle_fastqs:
    input:
        f"{outdir}/unshuffled/{{sample}}_R{{p}}.fq.gz"
        if config["seq_tech"] == "illumina"
        else f"{outdir}/unshuffled/{{sample}}.fq.gz",
    output:
        f"{outdir}/fastq/{{sample}}_R{{p}}.fq.gz"
        if config["seq_tech"] == "illumina"
        else f"{outdir}/fastq/{{sample}}.fq.gz",
    log:
        f"{outdir}/logs/seqkit/{{sample}}_R{{p}}.log"
        if config["seq_tech"] == "illumina"
        else f"{outdir}/logs/seqkit/{{sample}}.log",
    params:
        lambda wildcards: shuf_seed[wildcards.sample],
    threads: 3
    shell:
        """
        seqkit \\
        shuffle {input} \\
        -j {threads} \\
        -s {params} \\
        -o {output} 2> {log}
        """
