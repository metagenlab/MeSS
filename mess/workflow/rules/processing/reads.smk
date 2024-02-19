if SEQ_TECH == "illumina":

    rule convert_sam_to_bam:
        input:
            sam_out,
        output:
            temp(os.path.join(dir.out.bam, "{sample}", "{fasta}.bam")),
        benchmark:
            os.path.join(dir.out.bench, "samtools", "sam2bam", "{sample}-{fasta}.txt")
        log:
            os.path.join(dir.out.logs, "samtools", "sam2bam", "{sample}-{fasta}.log"),
        resources:
            mem_mb=config.resources.norm.mem,
            mem=str(config.resources.norm.mem) + "MB",
            time=config.resources.norm.time,
        threads: config.resources.norm.cpu
        conda:
            os.path.join(dir.env, "bioconvert.yml")
        shell:
            """
            samtools view -@ {threads} -bS {input} | \
            samtools sort -@ {threads} > {output} 2> {log}
            """

    fastq_gz = [
        temp(os.path.join(dir.out.short, "{sample}", "{fasta}1.fq.gz")),
        temp(os.path.join(dir.out.short, "{sample}", "{fasta}2.fq.gz")),
    ]

    rule compress_fastq:
        input:
            fastq_out,
        output:
            fastq_gz if PAIRED else fastq_gz[0],
        params:
            prefix=fq_prefix,
        benchmark:
            os.path.join(dir.out.bench, "pigz", "{sample}-{fasta}.txt")
        resources:
            mem_mb=config.resources.norm.mem,
            mem=str(config.resources.norm.mem) + "MB",
            time=config.resources.med.time,
        threads: config.resources.norm.cpu
        shell:
            """
            pigz -p {threads} {params.prefix}*.fq
            """


rule merge_bams:
    input:
        lambda wildcards: list_cat(wildcards, "bam"),
    output:
        os.path.join(dir.out.bam, "{sample}.bam"),
    benchmark:
        os.path.join(dir.out.bench, "samtools", "merge", "{sample}.txt")
    log:
        os.path.join(dir.out.logs, "samtools", "merge", "{sample}.log"),
    resources:
        mem_mb=config.resources.norm.mem,
        mem=str(config.resources.norm.mem) + "MB",
        time=config.resources.norm.time,
    threads: config.resources.norm.cpu
    conda:
        os.path.join(dir.env, "bioconvert.yml")
    shell:
        """
        samtools merge -@ {threads} -o {output} {input} 2> {log}
        """


rule index_bam:
    input:
        os.path.join(dir.out.bam, "{sample}.bam"),
    output:
        os.path.join(dir.out.bam, "{sample}.bam.bai"),
    benchmark:
        os.path.join(dir.out.bench, "samtools", "index", "{sample}.txt")
    resources:
        mem_mb=config.resources.norm.mem,
        mem=str(config.resources.norm.mem) + "MB",
        time=config.resources.norm.time,
    threads: config.resources.norm.cpu
    conda:
        os.path.join(dir.env, "bioconvert.yml")
    shell:
        """
        samtools index -@ {threads} {input}
        """


sample_fastq_out = []
if SKIP_SHUFFLE:
    if SEQ_TECH == "illumina":
        sample_fastq_out = os.path.join(dir.out.fastq, "{sample}_R{p}.fq.gz")
    else:
        sample_fastq_out = os.path.join(dir.out.fastq, "{sample}.fq.gz")
else:
    if SEQ_TECH == "illumina":
        sample_fastq_out = temp(os.path.join(dir.out.cat, "{sample}_R{p}.fq.gz"))
    else:
        sample_fastq_out = temp(os.path.join(dir.out.cat, "{sample}.fq.gz"))


rule cat_fastqs:
    input:
        lambda wildcards: list_cat(wildcards, "fastq"),
    output:
        sample_fastq_out,
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.norm.time,
    shell:
        """
        cat {input} > {output}
        """


if not SKIP_SHUFFLE:

    rule shuffle_fastqs:
        input:
            os.path.join(dir.out.cat, "{sample}_R{p}.fq.gz")
            if SEQ_TECH == "illumina"
            else os.path.join(dir.out.cat, "{sample}.fq.gz"),
        output:
            temp(os.path.join(dir.out.shuffle, "{sample}_R{p}.fq.gz"))
            if SEQ_TECH == "illumina"
            else temp(os.path.join(dir.out.shuffle, "{sample}.fq.gz")),
        params:
            lambda wildcards: SHUFFLE[wildcards.sample],
        benchmark:
            os.path.join(
                dir.out.bench, "seqkit", "shuffle", "{sample}_R{p}.txt"
            ) if SEQ_TECH == "illumina" else os.path.join(
                dir.out.bench, "seqkit", "shuffle", "{sample}.txt"
            )
        log:
            os.path.join(dir.out.logs, "seqkit", "shuffle", "{sample}_R{p}.log")
            if SEQ_TECH == "illumina"
            else os.path.join(dir.out.logs, "seqkit", "shuffle", "{sample}.log"),
        resources:
            mem_mb=config.resources.med.mem,
            mem=str(config.resources.med.mem) + "MB",
            time=config.resources.norm.time,
        threads: config.resources.norm.cpu
        conda:
            os.path.join(dir.env, "seqkit.yml")
        shell:
            """
            seqkit \\
            shuffle {input} \\
            -j {threads} \\
            -s {params} \\
            -o {output} 2> {log}
            """

    rule anonymize_reads:
        input:
            os.path.join(dir.out.shuffle, "{sample}_R{p}.fq.gz")
            if SEQ_TECH == "illumina"
            else os.path.join(dir.out.shuffle, "{sample}.fq.gz"),
        output:
            os.path.join(dir.out.fastq, "{sample}_R{p}.fq.gz")
            if SEQ_TECH == "illumina"
            else os.path.join(dir.out.fastq, "{sample}.fq.gz"),
        conda:
            os.path.join(dir.env, "seqkit.yml")
        benchmark:
            os.path.join(
                dir.out.bench, "seqkit", "anonymize", "{sample}_R{p}.txt"
            ) if SEQ_TECH == "illumina" else os.path.join(
                dir.out.bench, "seqkit", "anonymize", "{sample}.txt"
            )
        params:
            lambda wildcards: f"{wildcards.sample}_{{nr}}/{wildcards.p}"
            if SEQ_TECH == "illumina"
            else lambda wildcards: f"{wildcards.sample}_{{nr}}",
        log:
            os.path.join(dir.out.logs, "seqkit", "replace", "{sample}_R{p}.log")
            if SEQ_TECH == "illumina"
            else os.path.join(dir.out.logs, "seqkit", "replace", "{sample}.log"),
            os.path.join(dir.out.logs, "anonkeys", "{sample}_R{p}.tsv"),
        resources:
            mem_mb=config.resources.norm.mem,
            mem=str(config.resources.norm.mem) + "MB",
            time=config.resources.norm.time,
        shell:
            """
            zcat {input} | seqkit replace \\
            -p .+ -r "{params}" -o {output} 2> {log[0]}
            paste -d '\t' <(seqkit seq -n {output}) <(seqkit seq -n {input}) > {log[1]} 
            """
