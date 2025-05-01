fastq_dir = dir.out.long
contig = "{contig}"
if CIRCULAR:
    contig = "{contig}_{n}"
sam_in = os.path.join(dir.out.bam, "{sample}", "{fasta}", contig + ".sam")
sam_in_ef = os.path.join(dir.out.ef, "{sample}", "{fasta}", contig + ".sam")

if SEQ_TECH == "illumina":
    fastq_dir = dir.out.short
    sam_in = os.path.join(fastq_dir, "{sample}", "{fasta}", contig + ".sam")
    sam_in_ef = (
        os.path.join(fastq_dir, "{sample}", "{fasta}", contig + "_errFree.sam"),
    )

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
            os.path.join(
                dir.out.logs, "sam2bam", "{sample}", "{fasta}", contig + ".ccs.log"
            ),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        threads: config.resources.norm.cpu
        conda:
            os.path.join(dir.conda, "samtools.yml")
        container:
            containers.samtools
        shell:
            """
            samtools view -@ {threads} -Sb {input} | \\
            samtools sort -@ {threads} -o {output} 2> {log}
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
if BAM or TAX:

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


rule convert_sam_to_bam:
    input:
        sam_in,
    output:
        temp(os.path.join(dir.out.bam, "{sample}", "{fasta}", contig + ".bam")),
    log:
        os.path.join(dir.out.logs, "sam2bam", "{sample}", "{fasta}", contig + ".log"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.norm.cpu
    conda:
        os.path.join(dir.conda, "samtools.yml")
    container:
        containers.samtools
    shell:
        """
        samtools view -@ {threads} -Sb {input} | \\
        samtools sort -@ {threads} -o {output} 2> {log}
        """


rule merge_bams:
    input:
        lambda wildcards: aggregate(wildcards, dir.out.bam, "bam"),
    output:
        os.path.join(dir.out.bam, "{sample}.bam"),
    log:
        os.path.join(
            dir.out.logs,
            "merge",
            "{sample}.log",
        ),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.norm.cpu
    conda:
        os.path.join(dir.conda, "samtools.yml")
    container:
        containers.samtools
    shell:
        """
        samtools merge -@ {threads} -o {output} {input} 2> {log}
        """


rule index_bams:
    input:
        os.path.join(dir.out.bam, "{sample}.bam"),
    output:
        os.path.join(dir.out.bam, "{sample}.bam.bai"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.norm.time,
    threads: config.resources.norm.cpu
    conda:
        os.path.join(dir.conda, "samtools.yml")
    container:
        containers.samtools
    shell:
        """
        samtools index -@ {threads} {input}
        """


rule get_bam_coverage:
    input:
        os.path.join(dir.out.bam, "{sample}.bam"),
    output:
        temp(os.path.join(dir.out.bam, "{sample}.txt")),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    conda:
        os.path.join(dir.conda, "samtools.yml")
    container:
        containers.samtools
    shell:
        """
        samtools coverage {input} > {output}
        """


rule get_contig_taxonomy:
    input:
        tax=os.path.join(dir.out.base, "coverages.tsv"),
        cov=os.path.join(dir.out.processing, "cov.tsv"),
    output:
        os.path.join(dir.out.processing, "tax.tsv"),
    run:
        taxdf = pd.read_csv(input.tax, sep="\t")
        covdf = pd.read_csv(input.cov, sep="\t")
        same_cols = list(np.intersect1d(taxdf.columns, covdf.columns))
        covdf.merge(taxdf, on=same_cols, how="left").to_csv(
            output[0], sep="\t", index=False
        )


rule get_tax_profile:
    input:
        cov=os.path.join(dir.out.bam, "{sample}.txt"),
        tax=os.path.join(dir.out.processing, "tax.tsv"),
    output:
        counts=os.path.join(dir.out.tax, "{sample}.tsv"),
        seq_abundance=os.path.join(dir.out.tax, "{sample}_seq.tsv"),
        tax_abundance=os.path.join(dir.out.tax, "{sample}_tax.tsv"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    params:
        circular=CIRCULAR,
    run:
        tax_df = pd.read_csv(input.tax, sep="\t")
        tax_df = tax_df[tax_df.samplename == wildcards.sample]
        if params.circular:
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
            df[["tax_id", "abundance"]].to_csv(out, sep="\t", header=False, index=False)


rule get_cami_profile:
    input:
        tsv=os.path.join(dir.out.tax, "{sample}_{abundance}.tsv"),
        dmp=os.path.join(TAXONKIT, "names.dmp"),
    output:
        os.path.join(dir.out.tax, "{sample}_{abundance}_CAMI.txt"),
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


rule cat_fastqs:
    input:
        lambda wildcards: aggregate(wildcards, fastq_dir, "fq.gz"),
    output:
        sample_fastq_out,
    params:
        dir=os.path.join(fastq_dir, "{sample}"),
        name="*{p}.fq.gz" if PAIRED else "*.fq.gz",
        head=lambda wildcards, input: list(input)[:3],
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.norm.time,
    message:
        "Concatenating {wildcards.sample} reads : {params.head} ... "
    shell:
        """
        find {params.dir} -name "{params.name}" | xargs cat > {output}
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


if ERRFREE:

    rule convert_sam_to_bam_ef:
        input:
            sam_in_ef,
        output:
            temp(os.path.join(dir.out.ef, "{sample}", "{fasta}", contig + ".bam")),
        log:
            os.path.join(
                dir.out.logs, "sam2bam", "{sample}", "{fasta}", contig + "_ef.log"
            ),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        threads: config.resources.norm.cpu
        conda:
            os.path.join(dir.conda, "samtools.yml")
        container:
            containers.samtools
        shell:
            """
            samtools view -@ {threads} -Sb {input} | \\
            samtools sort -@ {threads} -o {output} 2> {log}
            """

    rule merge_bams_ef:
        input:
            lambda wildcards: aggregate(wildcards, dir.out.ef, "bam"),
        output:
            os.path.join(dir.out.ef, "{sample}.bam"),
        log:
            os.path.join(
                dir.out.logs,
                "samtools",
                "merge",
                "{sample}_ef.log",
            ),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        threads: config.resources.norm.cpu
        conda:
            os.path.join(dir.conda, "samtools.yml")
        container:
            containers.samtools
        shell:
            """
            samtools merge -@ {threads} -o {output} {input} 2> {log}
            """

    rule index_bams_ef:
        input:
            os.path.join(dir.out.ef, "{sample}.bam"),
        output:
            os.path.join(dir.out.ef, "{sample}.bam.bai"),
        benchmark:
            os.path.join(dir.out.bench, "samtools", "index", "{sample}_ef.txt")
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.norm.time,
        threads: config.resources.norm.cpu
        conda:
            os.path.join(dir.conda, "samtools.yml")
        container:
            containers.samtools
        shell:
            """
            samtools index -@ {threads} {input}
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
