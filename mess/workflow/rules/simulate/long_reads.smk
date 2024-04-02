fa = os.path.join(dir.out.base, "split", "{fasta}", "{contig}.fa")


rule pbsim3:
    input:
        fa=fa,
        flag=os.path.join(dir.out.base, "split.tsv"),
        df=os.path.join(dir.out.base, "cov.tsv"),
    output:
        temp(os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}_0001.sam"))
        if PASSES > 1
        else temp(
            os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}_0001.fastq")
        ),
        temp(os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}_0001.maf")),
        temp(os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}_0001.ref")),
    params:
        model=MODEL,
        ratio=RATIO,
        meanlen=MEAN_LEN,
        lensd=SD_LEN,
        minlen=MIN_LEN,
        maxlen=MAX_LEN,
        passes=PASSES,
        accuracy=ACCURACY,
        cov=lambda wildcards, input: get_value(input.df, wildcards, "cov_sim"),
        seed=lambda wildcards, input: int(get_value(input.df, wildcards, "seed")),
        prefix=os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}"),
    benchmark:
        os.path.join(dir.out.bench, "pbsim3", "{sample}", "{fasta}", "{contig}.txt")
    log:
        os.path.join(dir.out.logs, "pbsim3", "{sample}", "{fasta}", "{contig}.log"),
    resources:
        mem_mb=config.resources.norm.mem,
        mem=str(config.resources.norm.mem) + "MB",
        time=config.resources.norm.time,
    conda:
        os.path.join(dir.env, "pbsim3.yml")
    shell:
        """
        pbsim --strategy wgs --method qshmm \\
        --id-prefix {wildcards.contig} \\
        --difference-ratio {params.ratio} \\
        --length-min {params.minlen} \\
        --length-max {params.maxlen} \\
        --length-mean {params.meanlen} \\
        --length-sd {params.lensd} \\
        --prefix {params.prefix} \\
        --seed {params.seed} \\
        --qshmm ${{CONDA_PREFIX}}/data/{params.model}.model \\
        --pass-num {params.passes} \\
        --accuracy-mean {params.accuracy} \\
        --depth {params.cov} --genome {input.fa} &> {log}
        """


tmpdir = os.path.join(OUTPUT, "tmp")

if PASSES > 1:

    rule ccs_sam_to_bam:
        input:
            os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}_0001.sam"),
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
        params:
            passes=PASSES,
            accuracy=ACCURACY,
        resources:
            mem_mb=config.resources.med.mem,
            mem=str(config.resources.med.mem) + "MB",
            time=config.resources.med.time,
            tmpdir=tmpdir,
        benchmark:
            os.path.join(dir.out.bench, "ccs", "{sample}", "{fasta}", "{contig}.txt")
        log:
            os.path.join(dir.out.logs, "ccs", "{sample}", "{fasta}", "{contig}.log"),
        conda:
            os.path.join(dir.env, "pbccs.yml")
        shell:
            """
            MIMALLOC_PAGE_RESET=0 MIMALLOC_LARGE_OS_PAGES=1 \\
            ccs -j {threads} --min-passes {params.passes} \\
            --min-rq {params.accuracy} {input} {output.fq} --report-file {log}
            """


if BAM:

    rule add_reference_name:
        input:
            maf=os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}_0001.maf"),
            fa=os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}_0001.ref"),
        output:
            temp(os.path.join(dir.out.bam, "{sample}", "{fasta}", "{contig}.maf")),
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
            mem_mb=config.resources.norm.mem,
            mem=str(config.resources.norm.mem) + "MB",
            time=config.resources.norm.time,
        threads: config.resources.sml.cpu
        conda:
            os.path.join(dir.env, "bioconvert.yml")
        shell:
            """
            bioconvert {input} {output} 2>> {log}
            """


if PASSES == 1:

    rule rename_single_pass_fastq:
        input:
            os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}_0001.fastq"),
        output:
            temp(os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}.fq")),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        threads: config.resources.sml.cpu
        shell:
            """
            mv {input} {output}
            """
