fastq_dir = dir.out.long
sam_in = os.path.join(dir.out.bam, "{sample}", "{fasta}", "{contig}.sam")
if SEQ_TECH == "illumina":
    fastq_dir = dir.out.short
    sam_in = os.path.join(fastq_dir, "{sample}", "{fasta}", "{contig}.sam")

fastq = os.path.join(fastq_dir, "{sample}", "{fasta}", "{contig}.fq")
fastq_gz = temp(os.path.join(fastq_dir, "{sample}", "{fasta}", "{contig}.fq.gz"))
if PAIRED:
    fastq = os.path.join(fastq_dir, "{sample}", "{fasta}", "{contig}{p}.fq")
    fastq_gz = temp(os.path.join(fastq_dir, "{sample}", "{fasta}", "{contig}{p}.fq.gz"))

# hifi reads post processing
if PASSES > 1:

    rule ccs_sam_to_bam:
        input:
            os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}.sam"),
        output:
            temp(
                os.path.join(
                    dir.out.base,
                    "ccs",
                    "{sample}",
                    "{fasta}",
                    "{contig}.ccs.bam",
                )
            ),
        benchmark:
            os.path.join(
                dir.out.bench,
                "samtools",
                "sam2bam",
                "{sample}",
                "{fasta}_{contig}.txt",
            )
        log:
            os.path.join(dir.out.logs, "ccs", "{sample}", "{fasta}", "{contig}.log"),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.norm.time,
        threads: config.resources.norm.cpu
        conda:
            os.path.join(dir.env, "bioconvert.yml")
        shell:
            """
            samtools view -@ {threads} -bS {input} | \
            samtools sort -@ {threads} > {output} 2> {log}
            """

    rule ccs_bam_to_fastq:
        input:
            os.path.join(dir.out.base, "ccs", "{sample}", "{fasta}", "{contig}.ccs.bam"),
        output:
            fq=temp(
                os.path.join(
                    dir.out.long,
                    "{sample}",
                    "{fasta}",
                    "{contig}.fq.gz",
                )
            ),
            json=temp(
                os.path.join(
                    dir.out.long,
                    "{sample}",
                    "{fasta}",
                    "{contig}.zmw_metrics.json.gz",
                )
            ),
        benchmark:
            os.path.join(dir.out.bench, "ccs", "{sample}", "{fasta}", "{contig}.txt")
        log:
            os.path.join(dir.out.logs, "ccs", "{sample}", "{fasta}", "{contig}.log"),
        params:
            passes=PASSES,
            accuracy=ACCURACY,
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.norm.time,
        threads: config.resources.norm.cpu
        conda:
            os.path.join(dir.env, "pbccs.yml")
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
            maf=os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}.maf"),
            fa=os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}.ref"),
        output:
            temp(os.path.join(dir.out.bam, "{sample}", "{fasta}", "{contig}.maf")),
        params:
            seqname=lambda wildcards, input: get_header(input.fa),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        conda:
            os.path.join(dir.env, "utils.yml")
        shell:
            """
            sed 's/ref/{params.seqname}/g' {input.maf} > {output}
            """

    rule convert_maf_to_sam:
        input:
            os.path.join(dir.out.bam, "{sample}", "{fasta}", "{contig}.maf"),
        output:
            temp(os.path.join(dir.out.bam, "{sample}", "{fasta}", "{contig}.sam")),
        benchmark:
            os.path.join(
                dir.out.logs,
                "bioconvert",
                "maf2sam",
                "{sample}",
                "{fasta}_{contig}.txt",
            )
        log:
            os.path.join(
                dir.out.logs,
                "bioconvert",
                "maf2sam",
                "{sample}",
                "{fasta}_{contig}.log",
            ),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        threads: config.resources.sml.cpu
        conda:
            os.path.join(dir.env, "bioconvert.yml")
        shell:
            """
            bioconvert {input} {output} 2> {log}
            """


rule convert_sam_to_bam:
    input:
        sam_in,
    output:
        temp(os.path.join(dir.out.bam, "{sample}", "{fasta}", "{contig}.bam")),
    log:
        os.path.join(
            dir.out.logs, "bioconvert", "sam2bam", "{sample}", "{fasta}{contig}.log"
        ),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    conda:
        os.path.join(dir.env, "bioconvert.yml")
    shell:
        """
        bioconvert {input} {output} -t {threads} 2> {log}
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
        os.path.join(dir.env, "bioconvert.yml")
    shell:
        """
        samtools merge -@ {threads} -o {output} {input} 2> {log}
        """


rule get_sample_bams:
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
        os.path.join(dir.env, "bioconvert.yml")
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
        os.path.join(dir.env, "bioconvert.yml")
    shell:
        """
        samtools sort -@ {threads} {input} -o {output}  2> {log}
        """


if BAM:

    rule download_taxdump:
        output:
            temp(os.path.join(TAXONKIT, "taxdump.tar.gz")),
        log:
            os.path.join(dir.out.logs, "curl", "taxdump.log"),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        conda:
            os.path.join(dir.env, "utils.yml")
        shell:
            """
            curl https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz \\
            -o {output} 2> {log}
            """

    rule decompress_taxdump:
        input:
            os.path.join(TAXONKIT, "taxdump.tar.gz"),
        output:
            os.path.join(TAXONKIT, "names.dmp"),
            os.path.join(TAXONKIT, "nodes.dmp"),
            os.path.join(TAXONKIT, "delnodes.dmp"),
            os.path.join(TAXONKIT, "merged.dmp"),
        params:
            dir=TAXONKIT,
        log:
            os.path.join(dir.out.logs, "tar", "taxdump.log"),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        log:
            os.path.join(dir.out.logs, "wget", "taxdump.log"),
        conda:
            os.path.join(dir.env, "utils.yml")
        shell:
            """
            tar -xzvf {input} -C {params.dir} &> {log}
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
        os.path.join(dir.env, "bioconvert.yml")
    shell:
        """
        samtools coverage {input} > {output}
        """


rule get_tax_profile:
    input:
        cov=os.path.join(dir.out.bam, "{sample}.txt"),
        tax=get_cov_table,
    output:
        temp(os.path.join(dir.out.tax, "{sample}.tsv")),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    run:
        tax_df = pd.read_csv(input.tax, sep="\t")
        cov_df = pd.read_csv(input.cov, sep="\t")
        cov_df.rename(columns={"#rname": "contig"}, inplace=True)
        merge_df = tax_df.merge(cov_df)
        df = merge_df.groupby("tax_id")["meandepth"].mean().reset_index()
        df["tax_abundance"] = df["meandepth"] / df["meandepth"].sum()
        df[["tax_id", "tax_abundance"]].to_csv(
            output[0], sep="\t", header=False, index=False
        )


rule tax_profile_to_biobox:
    input:
        tsv=os.path.join(dir.out.tax, "{sample}.tsv"),
        dmp=os.path.join(TAXONKIT, "names.dmp"),
    output:
        os.path.join(dir.out.tax, "{sample}_profile.txt"),
    log:
        os.path.join(dir.out.logs, "taxonkit", "profile2cami", "{sample}.log"),
    params:
        dir=TAXONKIT,
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    conda:
        os.path.join(dir.env, "taxonkit.yml")
    shell:
        """
        taxonkit \\
            profile2cami \\
            -j {threads} \\
            --data-dir {params.dir} \\
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
        os.path.join(dir.env, "bioconvert.yml")
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
        os.path.join(dir.env, "utils.yml")
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
            os.path.join(
                dir.out.bench, "seqkit", "shuffle", "{sample}_R{p}.txt"
            ) if PAIRED else os.path.join(
                dir.out.bench, "seqkit", "shuffle", "{sample}.txt"
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
            if PAIRED
            else os.path.join(dir.out.shuffle, "{sample}.fq.gz"),
        output:
            os.path.join(dir.out.fastq, "{sample}_R{p}.fq.gz")
            if PAIRED
            else os.path.join(dir.out.fastq, "{sample}.fq.gz"),
        conda:
            os.path.join(dir.env, "seqkit.yml")
        benchmark:
            os.path.join(
                dir.out.bench, "seqkit", "anonymize", "{sample}_R{p}.txt"
            ) if PAIRED else os.path.join(
                dir.out.bench, "seqkit", "anonymize", "{sample}.txt"
            )
        params:
            seqkit_replace,
        log:
            os.path.join(dir.out.logs, "seqkit", "replace", "{sample}_R{p}.log")
            if PAIRED
            else os.path.join(dir.out.logs, "seqkit", "replace", "{sample}.log"),
            os.path.join(dir.out.logs, "anonymized_reads", "{sample}_R{p}.tsv")
            if PAIRED
            else os.path.join(dir.out.logs, "anonymized_reads", "{sample}.tsv"),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.norm.time,
        shell:
            """
            seqkit seq {input} | seqkit replace \\
            -p .+ -r "{params}" -o {output} 2> {log[0]}
            paste -d '\t' <(seqkit seq -n {output}) <(seqkit seq -n {input}) > {log[1]} 
            """
