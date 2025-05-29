rule rename_fastas:
    input:
        fasta_input,
    output:
        temp(os.path.join(dir.out.processing, "{fasta}.fasta")),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    conda:
        os.path.join(dir.conda, "seqkit.yml")
    container:
        containers.seqkit
    shell:
        """
        seqkit \\
        seq \\
        -i {input} \\
        > {output}
        """


if AMPLICONS:

    if PRIMERS:
        primers = os.path.abspath(PRIMERS)
    else:
        primers = os.path.join(dir.out.processing, "primers.txt")

    rule primersearch:
        input:
            primers=primers,
            fasta=os.path.join(dir.out.processing, "{fasta}.fasta"),
        output:
            temp(os.path.join(dir.out.processing, "{fasta}.primersearch.txt")),
        params:
            mismatch=MISMATCH,
        log:
            os.path.join(dir.out.logs, "primersearch", "{fasta}.log"),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        conda:
            os.path.join(dir.conda, "emboss.yml")
        container:
            containers.emboss
        shell:
            """
            primersearch \\
            -seqall {input.fasta} \\
            -infile {input.primers} \\
            -outfile {output} \\
            -mismatchpercent {params.mismatch} \\
            2> {log}
            """

    amp_args = ""
    if CUT:
        amp_args += "--cut "
    if ORIENT:
        amp_args += "--orient "

    rule get_amplicons:
        input:
            search=os.path.join(dir.out.processing, "{fasta}.primersearch.txt"),
            fasta=os.path.join(dir.out.processing, "{fasta}.fasta"),
            primers=primers,
        output:
            temp(os.path.join(dir.out.processing, "{fasta}.amplicons.fasta")),
        params:
            minlen=AMP_MINLEN,
            maxlen=AMP_MAXLEN,
            args=amp_args,
            script=os.path.join(dir.scripts, "get_amplicons.py"),
        log:
            os.path.join(dir.out.logs, "primersearch", "{fasta}_summary.log"),
        shell:
            """
            python {params.script} \\
            --input {input.search} \\
            --fasta {input.fasta} \\
            --primers {input.primers} \\
            --minlen {params.minlen} \\
            --maxlen {params.maxlen} \\
            {params.args} \\
            --output {output} \\
            --log 2> {log}
            """


rule get_fasta_stats:
    input:
        list_fastas,
    output:
        os.path.join(dir.out.processing, "seqkit_stats.tsv"),
    params:
        os.path.join(dir.out.processing, "*.fasta")
        if not AMPLICONS
        else os.path.join(dir.out.processing, "*.amplicons.fasta"),
    log:
        os.path.join(dir.out.logs, "seqkit", "stats.log"),
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
        stats -b -T -j {threads} \\
        {params} \\
        > {output} \\
        2> {log}
        """


checkpoint split_contigs:
    input:
        fa=list_fastas,
        cov=os.path.join(dir.out.base, "coverages.tsv"),
    output:
        tsv=os.path.join(dir.out.processing, "cov.tsv"),
        dir=directory(os.path.join(dir.out.processing, "split")),
    params:
        circular=CIRCULAR,
        rotate=ROTATE,
        amplicons=AMPLICONS,
        read_len=MEAN_LEN,
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    script:
        os.path.join(dir.scripts, "split_contigs.py")


if CIRCULAR:

    rule rotate_contigs:
        input:
            os.path.join(dir.out.processing, "split", "{fasta}_{contig}.fna"),
        output:
            os.path.join(
                dir.out.processing, "rotate", "{sample}", "{fasta}_{contig}_{n}.fna"
            ),
        params:
            lambda wildcards: get_value("random_start", wildcards),
        log:
            os.path.join(
                dir.out.logs,
                "seqkit",
                "restart",
                "{sample}",
                "{fasta}_{contig}_{n}.log",
            ),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        threads: config.resources.sml.cpu
        conda:
            os.path.join(dir.conda, "seqkit.yml")
        container:
            containers.seqkit
        shell:
            """
            seqkit restart \\
            -i {params} \\
            {input} | \\
            seqkit replace -p .+ -r {wildcards.contig}_{wildcards.n} \\
            > {output}
            """
