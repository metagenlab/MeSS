rule pbsim3:
    input:
        fasta=os.path.join(dir.out.fasta, "{fasta}.renamed"),
        flag=os.path.join(dir.out.fasta, "split", "split.tsv"),
        fa=os.path.join(dir.out.fasta, "split", "{fasta}_{contig}.fa"),
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
    shell:
        """
        pbsim --strategy wgs --method qshmm \\
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
                    "{fasta}_{contig}.ccs.bam",
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
        shell:
            "samtools view -bS {input} | samtools sort > {output} 2> {log}"

    rule ccs_bam_to_fastq:
        input:
            os.path.join(dir.out.base, "ccs", "{sample}", "{fasta}", "{contig}.ccs.bam"),
        output:
            fq=temp(
                os.path.join(
                    dir.out.base,
                    "ccs",
                    "{sample}",
                    "{fasta}_{contig}_0001.fq.gz",
                )
            ),
            json=temp(
                os.path.join(
                    dir.out.base,
                    "ccs",
                    "{sample}",
                    "{fasta}_{contig}_0001.zmw_metrics.json.gz",
                )
            ),
        params:
            passes=PASSES,
            accuracy=ACCURACY,
        resources:
            tmpdir=tmpdir,
        benchmark:
            os.path.join(dir.out.bench, "ccs", "{sample}", "{fasta}", "{contig}.txt")
        log:
            os.path.join(dir.out.logs, "ccs", "{sample}", "{fasta}", "{contig}.log"),
        shell:
            """
            MIMALLOC_PAGE_RESET=0 MIMALLOC_LARGE_OS_PAGES=1 \\
            ccs --min-passes {params.passes} \\
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
        benchmark:
            os.path.join(dir.out.bench, "sed", "{sample}", "{fasta}", "{contig}.txt")
        shell:
            """
            sed 's/ref/{params.seqname}/g' {input.maf} > {output}
            """

    rule convert_maf_to_bam:
        input:
            os.path.join(dir.out.bam, "{sample}", "{fasta}", "{contig}.maf"),
        output:
            sam=temp(os.path.join(dir.out.bam, "{sample}", "{fasta}", "{contig}.sam")),
            bam=temp(os.path.join(dir.out.bam, "{sample}", "{fasta}", "{contig}.bam")),
        benchmark:
            os.path.join(
                dir.out.logs,
                "bioconvert",
                "maf2bam",
                "{sample}",
                "{fasta}_{contig}.txt",
            )
        log:
            os.path.join(
                dir.out.logs,
                "bioconvert",
                "maf2bam",
                "{sample}",
                "{fasta}_{contig}.log",
            ),
        threads: config.resources.med.cpu
        shell:
            """
            bioconvert {input} {output.sam} 2>> {log}
            samtools view -@ {threads} -bS {output.sam} | \\
            samtools sort -@ threads > {output.bam} 2>> {log}
            """

    rule merge_alignments_bam:
        input:
            lambda wildcards: pbsim3_expand(wildcards, dir.out.bam, "bam"),
        output:
            temp(os.path.join(dir.out.bam, "{sample}", "{fasta}.bam")),
        benchmark:
            os.path.join(dir.out.bench, "samtools", "merge", "{sample}", "{fasta}.txt")
        log:
            os.path.join(dir.out.logs, "samtools", "merge", "{sample}", "{fasta}.log"),
        threads: config.resources.med.cpu
        shell:
            "samtools merge -@ {threads} -o {output} {input} &>> {log}"


if PASSES == 1:

    rule compress_single_pass_fastq:
        input:
            os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}_0001.fastq"),
        output:
            temp(os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}.fq.gz")),
        params:
            os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}.fq"),
        threads: config.resources.med.cpu
        benchmark:
            os.path.join(dir.out.bench, "pigz", "{sample}", "{fasta}", "{contig}.txt")
        shell:
            """
            mv {input} {params}
            pigz -p {threads} {params}
            """


rule cat_contig_reads:
    input:
        lambda wildcards: pbsim3_expand(wildcards, dir.out.long, "fq.gz"),
    output:
        temp(os.path.join(dir.out.long, "{sample}", "{fasta}.fq.gz")),
    shell:
        "cat {input} > {output}"
