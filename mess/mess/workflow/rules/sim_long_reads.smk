multipass = False
if config["passes"] > 1:
    multipass = True
if config["seq_tech"] == "pacbio":
    model = "QSHMM-RSII"
    ratio = "6:55:39"
    if config["profile"] == "sequel":
        ratio = "22:45:33"
elif config["seq_tech"] == "nanopore":
    ratio = "39:24:36"
    if config["profile"] == "R10.4":
        model = "QSHMM-ONT-HQ"
    elif config["profile"] == "R10.3":
        model = "QSHMM-ONT"


wildcard_constraints:
    contig="[0-9]+",
    sample="[^/]+",
    fasta="[^/-]+",
    chunk="[0-9]+",


rule chunk_fasta:
    input:
        f"{outdir}/fasta/{{fasta}}.renamed",
    output:
        temp(f"{outdir}/fasta/{{fasta}}-{{chunk}}.fa"),
    shell:
        "cp {input} {output}"


rule pbsim3:
    input:
        fa=f"{outdir}/fasta/{{fasta}}-{{chunk}}.fa",
        df=f"{outdir}/cov.tsv",
    output:
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.sam")
        if multipass
        else temp(f"{outdir}/fastq/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.fastq"),
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.maf"),
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.ref"),
    params:
        model=model,
        ratio=ratio,
        meanlen=config["read_len"],
        lensd=config["read_sd"],
        minlen=config["min_len"],
        maxlen=config["max_len"],
        passes=config["passes"],
        accuracy=config["accuracy"],
        cov=lambda wildcards, input: get_value(input.df, wildcards, "cov_sim"),
        seed=lambda wildcards, input: int(get_value(input.df, wildcards, "seed")),
        prefix=f"{outdir}/fastq/{{sample}}/{{fasta}}-{{chunk}}",
    log:
        f"{outdir}/logs/pbsim3/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.log",
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


tmpdir = os.path.join(outdir, "tmp")

if multipass:

    rule ccs_sam_to_bam:
        input:
            f"{outdir}/fastq/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.sam",
        output:
            temp(f"{outdir}/ccs/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.ccs.bam"),
        log:
            f"{outdir}/logs/ccs/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.bam.log",
        shell:
            "samtools view -bS {input} | samtools sort > {output} 2> {log}"

    rule ccs_bam_to_fastq:
        input:
            f"{outdir}/ccs/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.ccs.bam",
        output:
            fq=temp(f"{outdir}/fastq/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.fq.gz"),
            json=temp(
                f"{outdir}/fastq/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.zmw_metrics.json.gz"
            ),
        log:
            f"{outdir}/logs/ccs/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.log",
        params:
            passes=config["passes"],
            accuracy=config["accuracy"],
        resources:
            tmpdir=tmpdir,
        shell:
            """
            MIMALLOC_PAGE_RESET=0 MIMALLOC_LARGE_OS_PAGES=1 \\
            ccs --min-passes {params.passes} \\
            --min-rq {params.accuracy} {input} {output.fq} --report-file {log}
            """


if config["bam"]:

    rule add_reference_name:
        input:
            maf=f"{outdir}/fastq/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.maf",
            fa=f"{outdir}/fastq/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.ref",
        output:
            temp(f"{outdir}/bam/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.maf"),
        params:
            seqname=lambda wildcards, input: get_header(input.fa),
        shell:
            """
            sed 's/ref/{params.seqname}/g' {input.maf} > {output}
            """

    rule convert_maf_to_bam:
        input:
            f"{outdir}/bam/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.maf",
        output:
            sam=temp(f"{outdir}/bam/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.sam"),
            bam=temp(f"{outdir}/bam/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.bam"),
        log:
            f"{outdir}/logs/bam/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.log",
        threads: 2
        shell:
            """
            bioconvert {input} {output.sam} 2>> {log}
            samtools view -@ {threads} -bS {output.sam} | \\
            samtools sort -@ threads > {output.bam} 2>> {log}
            """

    rule merge_alignments_bam:
        input:
            lambda wildcards: pbsim3_expand(wildcards, "bam", "bam", "both"),
        output:
            temp(f"{outdir}/bam/{{sample}}/{{fasta}}.bam"),
        log:
            f"{outdir}/logs/bam/{{sample}}/{{fasta}}.log",
        threads: 2
        shell:
            "samtools merge -@ {threads} -o {output} {input} &>> {log}"


if not multipass:

    rule compress_single_pass_fastq:
        input:
            f"{outdir}/fastq/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.fastq",
        output:
            temp(f"{outdir}/fastq/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.fq.gz"),
        params:
            f"{outdir}/fastq/{{sample}}/{{fasta}}-{{chunk}}_{{contig}}.fq",
        threads: 2
        shell:
            """
            mv {input} {params}
            pigz -p {threads} {params}
            """


rule cat_fastq_per_sample:
    input:
        lambda wildcards: pbsim3_expand(wildcards, "fastq", "fq.gz", "both"),
    output:
        temp(f"{outdir}/fastq/{{sample}}/{{fasta}}.fq.gz"),
    shell:
        "cat {input} > {output}"
