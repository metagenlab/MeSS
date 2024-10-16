fastq_dir = dir.out.long
contig = "{contig}"
if ROTATE > 1:
    contig = "{contig}_{n}"
sam_in = os.path.join(dir.out.bam, "{sample}", "{fasta}", contig + ".sam")
if SEQ_TECH == "illumina":
    fastq_dir = dir.out.short
    sam_in = os.path.join(fastq_dir, "{sample}", "{fasta}", contig + ".fixed")

fastq = os.path.join(fastq_dir, "{sample}", "{fasta}", "{contig}.fq")
fastq_gz = temp(os.path.join(fastq_dir, "{sample}", "{fasta}", "{contig}.fq.gz"))
if PAIRED:
    fastq = os.path.join(fastq_dir, "{sample}", "{fasta}", "{contig}{p}.fq")
    fastq_gz = temp(os.path.join(fastq_dir, "{sample}", "{fasta}", "{contig}{p}.fq.gz"))

# hifi reads post processing
if PASSES > 1:

    rule ccs_sam_to_bam:
        input:
            os.path.join(dir.out.long, "{sample}", "{fasta}", contig + ".sam"),
        output:
            temp(
                os.path.join(
                    dir.out.base,
                    "ccs",
                    "{sample}",
                    "{fasta}",
                    contig + ".ccs.bam",
                )
            ),
        log:
            os.path.join(dir.out.logs, "ccs", "{sample}", "{fasta}", contig + ".log"),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.norm.time,
        threads: config.resources.norm.cpu
        conda:
            os.path.join(dir.conda, "bioconvert.yml")
        container:
            containers.bioconvert
        shell:
            """
            samtools view -@ {threads} -bS {input} | \
            samtools sort -@ {threads} > {output} 2> {log}
            """

    rule ccs_bam_to_fastq:
        input:
            os.path.join(
                dir.out.base, "ccs", "{sample}", "{fasta}", contig + ".ccs.bam"
            ),
        output:
            fq=temp(
                os.path.join(
                    dir.out.long,
                    "{sample}",
                    "{fasta}",
                    contig + ".fq.gz",
                )
            ),
            json=temp(
                os.path.join(
                    dir.out.long,
                    "{sample}",
                    "{fasta}",
                    contig + ".zmw_metrics.json.gz",
                )
            ),
        log:
            os.path.join(dir.out.logs, "ccs", "{sample}", "{fasta}", contig + ".log"),
        params:
            passes=PASSES,
            accuracy=ACCURACY,
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.norm.time,
        threads: config.resources.norm.cpu
        conda:
            os.path.join(dir.conda, "pbccs.yml")
        container:
            containers.pbccs
        shell:
            """
            MIMALLOC_PAGE_RESET=0 MIMALLOC_LARGE_OS_PAGES=1 \\
            ccs -j {threads} --min-passes {params.passes} \\
            --min-rq {params.accuracy} {input} {output.fq} --report-file {log}
            """


# PBSIM3 bam pre-processing
if BAM:

    rule add_reference_name:
        input:
            maf=os.path.join(dir.out.long, "{sample}", "{fasta}", contig + ".maf"),
            fa=os.path.join(dir.out.long, "{sample}", "{fasta}", contig + ".ref"),
        output:
            temp(os.path.join(dir.out.bam, "{sample}", "{fasta}", contig + ".maf")),
        params:
            seqname=lambda wildcards, input: get_header(input.fa),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        shell:
            """
            sed 's/ref/{params.seqname}/g' {input.maf} > {output}
            """

    rule convert_maf_to_sam:
        input:
            os.path.join(dir.out.bam, "{sample}", "{fasta}", contig + ".maf"),
        output:
            temp(os.path.join(dir.out.bam, "{sample}", "{fasta}", contig + ".sam")),
        log:
            os.path.join(
                dir.out.logs,
                "bioconvert",
                "maf2sam",
                "{sample}",
                "{fasta}" + "_" + contig + ".log",
            ),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        threads: config.resources.sml.cpu
        conda:
            os.path.join(dir.conda, "bioconvert.yml")
        container:
            containers.bioconvert
        shell:
            """
            bioconvert {input} {output} 2> {log}
            """


rule fix_art_sam:
    """
    rule to replace SAM cigar string with read length + M
    Fixes truncated art_illumina SAM files with some genomes
    """
    input:
        os.path.join(fastq_dir, "{sample}", "{fasta}", contig + ".sam"),
    output:
        temp(os.path.join(fastq_dir, "{sample}", "{fasta}", contig + ".fixed")),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    params:
        maxlen=MEAN_LEN,
    shell:
        """
        awk 'BEGIN {{OFS="\t"}} {{ if ($1 ~ /^@/) {{ print $0 }} \\
        else {{ $6 = "{params.maxlen}M"; print $0 }} }}' \\
        {input} > {output}
        """


rule convert_sam_to_bam:
    input:
        sam_in,
    output:
        temp(os.path.join(dir.out.bam, "{sample}", "{fasta}", contig + ".bam")),
    log:
        os.path.join(
            dir.out.logs,
            "bioconvert",
            "sam2bam",
            "{sample}",
            "{fasta}" + contig + ".log",
        ),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    conda:
        os.path.join(dir.conda, "bioconvert.yml")
    container:
        containers.bioconvert
    shell:
        """
        bioconvert sam2bam {input} {output} -t {threads} 2> {log}
        """


rule merge_contig_bams:
    input:
        lambda wildcards: aggregate(wildcards, dir.out.bam, "contig", "bam"),
    output:
        temp(os.path.join(dir.out.bam, "{sample}", "{fasta}.bam")),
    benchmark:
        os.path.join(dir.out.bench, "samtools", "merge", "{sample}", "{fasta}.txt")
    log:
        os.path.join(dir.out.logs, "samtools", "merge", "{sample}", "{fasta}.log"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    conda:
        os.path.join(dir.conda, "bioconvert.yml")
    container:
        containers.bioconvert
    shell:
        """
        samtools merge -@ {threads} -o {output} {input} 2> {log}
        """


rule merge_sample_bams:
    input:
        lambda wildcards: aggregate(wildcards, dir.out.bam, "fasta", "bam"),
    output:
        temp(os.path.join(dir.out.bam, "{sample}.unsorted")),
    benchmark:
        os.path.join(dir.out.bench, "samtools", "merge", "{sample}.txt")
    log:
        os.path.join(dir.out.logs, "samtools", "merge", "{sample}.log"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.norm.cpu
    conda:
        os.path.join(dir.conda, "bioconvert.yml")
    container:
        containers.bioconvert
    shell:
        """
        samtools merge -@ {threads} -o {output} {input} 2> {log}
        """


rule sort_bams:
    input:
        os.path.join(dir.out.bam, "{sample}.unsorted"),
    output:
        os.path.join(dir.out.bam, "{sample}.bam"),
    benchmark:
        os.path.join(dir.out.bench, "samtools", "sort", "{sample}.txt")
    log:
        os.path.join(dir.out.logs, "samtools", "sort", "{sample}.log"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.norm.cpu
    conda:
        os.path.join(dir.conda, "bioconvert.yml")
    container:
        containers.bioconvert
    shell:
        """
        samtools sort -@ {threads} {input} -o {output} 2> {log}
        """


rule get_bam_coverage:
    input:
        os.path.join(dir.out.bam, "{sample}.bam"),
    output:
        temp(os.path.join(dir.out.bam, "{sample}.txt")),
    log:
        os.path.join(dir.out.logs, "samtools", "coverage", "{sample}.log"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    conda:
        os.path.join(dir.conda, "bioconvert.yml")
    container:
        containers.bioconvert
    shell:
        """
        samtools coverage {input} > {output} 2> {log}
        """


rule get_tax_profile:
    input:
        cov=os.path.join(dir.out.bam, "{sample}.txt"),
        tax=get_cov_table,
    output:
        counts=os.path.join(dir.out.tax, "{sample}.tsv"),
        seq_abundance=temp(os.path.join(dir.out.tax, "{sample}_seq.tsv")),
        tax_abundance=temp(os.path.join(dir.out.tax, "{sample}_tax.tsv")),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    params:
        rotate=ROTATE,
    run:
        tax_df = pd.read_csv(input.tax, sep="\t")
        tax_df = tax_df[tax_df.samplename == wildcards.sample]
        if params.rotate > 1:
            tax_df["contig"] = tax_df["contig"] + "_" + tax_df["n"].astype(str)
        cov_df = pd.read_csv(input.cov, sep="\t")
        cov_df.rename(columns={"#rname": "contig"}, inplace=True)
        merge_df = tax_df.merge(cov_df)
        merge_df[
            [
                "samplename",
                "fasta",
                "contig",
                "tax_id",
                "startpos",
                "endpos",
                "numreads",
                "covbases",
                "coverage",
                "cov_sim",
                "meandepth",
                "meanbaseq",
                "meanmapq",
            ]
        ].to_csv(output.counts, sep="\t", index=False)
        for col in ["numreads", "meandepth"]:
            if col == "numreads":
                out = output.seq_abundance
                df = merge_df.groupby("tax_id")[col].sum().reset_index()
            elif col == "meandepth":
                out = output.tax_abundance
                df = merge_df.groupby("tax_id")[col].mean().reset_index()
            df["abundance"] = df[col] / df[col].sum()
            df[["tax_id", "abundance"]].to_csv(
                out, sep="\t", header=False, index=False
            )



rule tax_profile_to_biobox:
    input:
        tsv=os.path.join(dir.out.tax, "{sample}_{abundance}.tsv"),
        dmp=os.path.join(TAXONKIT, "names.dmp"),
    output:
        os.path.join(dir.out.tax, "{sample}_{abundance}.txt"),
    log:
        os.path.join(
            dir.out.logs, "taxonkit", "profile2cami", "{sample}_{abundance}.log"
        ),
    params:
        dir=TAXONKIT,
        ranks=RANKS,
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    conda:
        os.path.join(dir.conda, "taxonkit.yml")
    container:
        containers.taxonkit
    shell:
        """
        taxonkit \\
        profile2cami \\
        -j {threads} --data-dir {params.dir} \\
        -r {params.ranks} \\
        -s {wildcards.sample} {input.tsv} > {output}
        """


rule index_bams:
    input:
        os.path.join(dir.out.bam, "{sample}.bam"),
    output:
        os.path.join(dir.out.bam, "{sample}.bam.bai"),
    benchmark:
        os.path.join(dir.out.bench, "samtools", "index", "{sample}.txt")
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.norm.time,
    threads: config.resources.norm.cpu
    conda:
        os.path.join(dir.conda, "bioconvert.yml")
    container:
        containers.bioconvert
    shell:
        """
        samtools index -@ {threads} {input}
        """


rule compress_contig_fastqs:
    input:
        fastq,
    output:
        fastq_gz,
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    conda:
        os.path.join(dir.conda, "pigz.yml")
    container:
        containers.pigz
    shell:
        """
        pigz -p {threads} {input}
        """


rule cat_contig_fastqs:
    input:
        flag=get_cov_table,
        fq=lambda wildcards: aggregate(wildcards, fastq_dir, "contig", "fq.gz"),
    output:
        temp(os.path.join(fastq_dir, "{sample}", "{fasta}{p}.fq.gz"))
        if PAIRED
        else temp(os.path.join(fastq_dir, "{sample}", "{fasta}.fq.gz")),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    shell:
        """
        cat {input.fq} > {output}
        """


sample_fastq_out = []
if SKIP_SHUFFLE:
    if PAIRED:
        sample_fastq_out = os.path.join(dir.out.fastq, "{sample}_R{p}.fq.gz")
    else:
        sample_fastq_out = os.path.join(dir.out.fastq, "{sample}.fq.gz")
else:
    if PAIRED:
        sample_fastq_out = temp(os.path.join(dir.out.cat, "{sample}_R{p}.fq.gz"))
    else:
        sample_fastq_out = temp(os.path.join(dir.out.cat, "{sample}.fq.gz"))


rule cat_sample_fastqs:
    input:
        lambda wildcards: aggregate(wildcards, fastq_dir, "fasta", "fq.gz"),
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
            if PAIRED
            else os.path.join(dir.out.cat, "{sample}.fq.gz"),
        output:
            temp(os.path.join(dir.out.shuffle, "{sample}_R{p}.fq.gz"))
            if PAIRED
            else temp(os.path.join(dir.out.shuffle, "{sample}.fq.gz")),
        params:
            lambda wildcards: SHUFFLE[wildcards.sample],
        benchmark:
            (
                os.path.join(dir.out.bench, "seqkit", "shuffle", "{sample}_R{p}.txt")
                if PAIRED
                else os.path.join(dir.out.bench, "seqkit", "shuffle", "{sample}.txt")
            )
        log:
            os.path.join(dir.out.logs, "seqkit", "shuffle", "{sample}_R{p}.log")
            if PAIRED
            else os.path.join(dir.out.logs, "seqkit", "shuffle", "{sample}.log"),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        threads: config.resources.norm.cpu
        conda:
            os.path.join(dir.conda, "seqkit.yml")
        container:
            containers.seqkit
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
            if PAIRED
            else os.path.join(dir.out.shuffle, "{sample}.fq.gz"),
        output:
            os.path.join(dir.out.fastq, "{sample}_R{p}.fq.gz")
            if PAIRED
            else os.path.join(dir.out.fastq, "{sample}.fq.gz"),
        benchmark:
            (
                os.path.join(dir.out.bench, "seqkit", "anonymize", "{sample}_R{p}.txt")
                if PAIRED
                else os.path.join(dir.out.bench, "seqkit", "anonymize", "{sample}.txt")
            )
        log:
            os.path.join(dir.out.logs, "seqkit", "replace", "{sample}_R{p}.log")
            if PAIRED
            else os.path.join(dir.out.logs, "seqkit", "replace", "{sample}.log"),
            os.path.join(dir.out.logs, "anonymized_reads", "{sample}_R{p}.tsv")
            if PAIRED
            else os.path.join(dir.out.logs, "anonymized_reads", "{sample}.tsv"),
        params:
            seqkit_replace,
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.norm.time,
        conda:
            os.path.join(dir.conda, "seqkit.yml")
        container:
            containers.bioconvert
        shell:
            """
            seqkit seq {input} | seqkit replace \\
            -p .+ -r "{params}" -o {output} 2> {log[0]}
            paste -d '\t' <(seqkit seq -n {output}) <(seqkit seq -n {input}) > {log[1]} 
            """


rule cleanup_files:
    input:
        list_reads,
        os.path.join(dir.out.base, "replicates.tsv"),
        os.path.join(dir.out.base, "samples.tsv"),
    output:
        temp(os.path.join(dir.out.base, "cleanup.done")),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.norm.time,
    params:
        proc=dir.out.processing,
        base=dir.out.base,
    shell:
        """
        rm -rf {params[0]}
        touch {output}
        """
