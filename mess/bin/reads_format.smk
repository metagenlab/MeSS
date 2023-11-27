rule convert_maf_to_sam:
    input:
        f"{outdir}/fastq/{{sample}}/{{fasta}}_0001.maf",
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
        if not config["paired"]
        else f"{outdir}/fastq/{{sample}}/{{fasta}}.sam",
    output:
        temp(f"{outdir}/bam/{{sample}}/{{fasta}}.bam"),
    log:
        f"{outdir}/logs/bam/{{sample}}/{{fasta}}.log",
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
        f"{outdir}/fastq/{{sample}}_R{{p}}.fq.gz"
        if config["seq_tech"] == "illumina"
        else f"{outdir}/fastq/{{sample}}.fq.gz",
    threads: 3
    shell:
        """
        cat {input} | pigz -p {threads} > {output}
        """
