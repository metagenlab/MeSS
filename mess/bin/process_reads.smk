rule convert_maf_to_sam:
    input:
        f"{outdir}/fastq/{{sample}}/{{fasta}}.maf",
    output:
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}.sam"),
    log:
        f"{outdir}/logs/bam/{{sample}}/{{fasta}}.log",
    shell:
        """
        bioconvert {input} {output} 2>> {log}
        """


rule convert_sam_to_bam:
    input:
        f"{outdir}/fastq/{{sample}}/{{fasta}}1.sam"
        if (not config["paired"]) and (config["seq_tech"]=="illumina")
        else f"{outdir}/fastq/{{sample}}/{{fasta}}.sam",
    output:
        temp(f"{outdir}/bam/{{sample}}/{{fasta}}.bam"),
    log:
        f"{outdir}/logs/bam/{{sample}}/{{fasta}}.log",
    shell:
        """
        bioconvert {input} {output} 2>> {log}
        """


rule concat_bam:
    input:
        lambda wildcards: list_concat(wildcards, "bam"),
    output:
        temp(f"{outdir}/bam/{{sample}}.unsorted"),
    threads: 3
    log:
        f"{outdir}/logs/bam/{{sample}}.log",
    shell:
        """
        samtools merge -@ {threads} -o {output} {input} 2>> {log}
        """


rule sort_bam:
    input:
        f"{outdir}/bam/{{sample}}.unsorted",
    output:
        f"{outdir}/bam/{{sample}}.bam",
    threads: 3
    log:
        f"{outdir}/logs/bam/{{sample}}.log",
    shell:
        """
        samtools sort -@ {threads} {input} > {output} 2>> {log}
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


rule concat_fastq:
    input:
        lambda wildcards: list_concat(wildcards, "fastq"),
    output:
        f"{outdir}/fastq/{{sample}}_R{{p}}.fq.gz"
        if config["seq_tech"] == "illumina"
        else f"{outdir}/fastq/{{sample}}.fq.gz",
    threads: 3
    shell:
        """
        cat {input} | pigz -p {threads} > {output}
        """
