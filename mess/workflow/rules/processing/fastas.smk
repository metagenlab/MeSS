rule symlink_fastas:
    output:
        temp(os.path.join(dir.out.processing, "links.txt")),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    run:
        get_fasta_table(wildcards).reset_index()[
            ["fasta", "path"]
        ].drop_duplicates().to_csv(output[0], sep="\t", header=False, index=False)


rule rename_fastas:
    input:
        os.path.join(dir.out.processing, "links.txt"),
    output:
        temp(os.path.join(dir.out.processing, "fasta", "renamed.done")),
    params:
        os.path.abspath(os.path.join(dir.out.processing, "fasta", "$1")),
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
        awk -F'\t' '{{print $2, $1}}' {input} | \\
        xargs -n 2 -P {threads} bash -c 'seqkit seq -i "$0" -o "{params}".fasta'
        touch {output}
        """


if PRIMERSEARCH:

    if PRIMERS:
        primers = os.path.abspath(PRIMERS)
    else:
        primers = os.path.join(dir.out.processing, "primers.txt")

    rule primersearch:
        input:
            primers=primers,
            renamed=os.path.join(dir.out.processing, "fasta", "renamed.done"),
        output:
            temp(os.path.join(dir.out.processing, "primersearch", "primersearch.done")),
        params:
            mismatch=MISMATCH,
            fa=os.path.abspath(os.path.join(dir.out.processing, "fasta", "*.fasta")),
            out=os.path.abspath(
                os.path.join(
                    dir.out.processing, "primersearch", "$base.primersearch.out"
                )
            ),
        log:
            os.path.join(dir.out.logs, "primersearch.log"),
        resources:
            mem_mb=config.resources.norm.mem,
            mem=str(config.resources.norm.mem) + "MB",
            time=config.resources.norm.time,
        threads: config.resources.norm.cpu
        conda:
            os.path.join(dir.conda, "emboss.yml")
        container:
            containers.emboss
        shell:
            """
            realpath {params.fa} | xargs -P {threads} -I {{}} sh -c '
            base=$(basename {{}} .fasta)
            primersearch \\
            -seqall {{}} \\
            -infile {input.primers} \\
            -outfile {params.out} \\
            -mismatchpercent {params.mismatch} 2>> {log}'
            touch {output}
            """

    rule get_primersearch_links:
        input:
            os.path.join(dir.out.processing, "primersearch", "primersearch.done"),
        output:
            os.path.join(dir.out.processing, "primersearch", "links.txt"),
        params:
            txt=os.path.abspath(
                os.path.join(dir.out.processing, "primersearch", "*.primersearch.out")
            ),
            fa=os.path.abspath(os.path.join(dir.out.processing, "fasta", "$base.fasta")),
        resources:
            mem_mb=config.resources.norm.mem,
            mem=str(config.resources.norm.mem) + "MB",
            time=config.resources.norm.time,
        threads: config.resources.norm.cpu
        shell:
            """
            realpath {params.txt} | \\
            xargs -P {threads} -I {{}} sh -c '
            base=$(basename {{}} .primersearch.out)
            echo "{{}}\\t{params.fa}\\t$base" >> {output}'
            """

    amp_args = ""
    if CUT:
        amp_args += "--cut "
    if KEEP_ORIENT:
        amp_args += "--keep-orient "

    rule get_amplicons:
        input:
            flag=os.path.join(dir.out.processing, "primersearch", "primersearch.done"),
            primers=primers,
            files=os.path.join(dir.out.processing, "primersearch", "links.txt"),
        output:
            temp(os.path.join(dir.out.processing, "fasta", "amplicons.done")),
        params:
            minlen=AMP_MINLEN,
            maxlen=AMP_MAXLEN,
            args=amp_args,
            script=os.path.join(dir.scripts, "get_amplicons.py"),
            out=os.path.abspath(
                os.path.join(dir.out.processing, "primersearch", "$2.fasta"),
            ),
        log:
            os.path.join(dir.out.logs, "primersearch.log"),
        shell:
            """
            awk -F'\\t' '{{print $1, $2, $3}}' {input.files} | \\
            xargs -n 3 -P {threads} bash -c '
            {params.script} \\
            -i $0 \\
            -f $1 \\
            -p {input.primers} \\
            -m {params.minlen} \\
            -M {params.maxlen} \\
            {params.args} \\
            --output {params.out} \\
            --log 2> {log}'
            touch {output}
            """


# rule get_seq_stats:
#     input:
#         os.path.join(dir.out.processing, "links.txt"),
#     output:
#         os.path.join(dir.out.processing, "seq_stats.tsv"),
#     log:
#         os.path.join(dir.out.logs, "seqkit", "stats.log"),
#     resources:
#         mem_mb=config.resources.sml.mem,
#         mem=str(config.resources.sml.mem) + "MB",
#         time=config.resources.sml.time,
#     threads: config.resources.norm.cpu
#     conda:
#         os.path.join(dir.conda, "seqkit.yml")
#     container:
#         containers.seqkit
#     shell:
#         """
#         cat {input} | \\
#         xargs -P {threads} -l bash -c 'seqkit fx2tab -n -i -l -g "$1" 2>> {log} | \\
#         awk -v fname="$0" "BEGIN{OFS=\"\t\"} {print fname, \$0}"' \\
#         > {output}
#         """


# checkpoint split_contigs:
#     input:
#         fa=list_fastas,
#         cov=os.path.join(dir.out.base, "coverages.tsv"),
#     output:
#         tsv=os.path.join(dir.out.processing, "cov.tsv"),
#         dir=directory(os.path.join(dir.out.processing, "split")),
#     params:
#         circular=CIRCULAR,
#         rotate=ROTATE,
#         read_len=MEAN_LEN,
#     resources:
#         mem_mb=config.resources.sml.mem,
#         mem=str(config.resources.sml.mem) + "MB",
#         time=config.resources.sml.time,
#     script:
#         os.path.join(dir.scripts, "split_contigs.py")


rule get_fa2ids:
    input:
        os.path.join(dir.out.processing, "links.txt"),
    output:
        temp(os.path.join(dir.out.processing, "fa2ids.tsv")),
    params:
        os.path.join(dir.out.processing, "fasta", "*.fasta"),
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
        awk -F'\\t' '{{print $2, $1}}' {input} | \\
        xargs -n 2 -P {threads} bash -c '
            seqkit fx2tab -nl "$1" \\
            | awk -v f=$(basename "$1") "{{print f\\\"\\\\t\\\"\\$0}}"
        ' _ > {output}
        """


if CIRCULAR:

    rule circularise_sequences:
        input:
            flag=os.path.join(dir.out.processing, "fasta", "renamed.done"),
            fa2id=os.path.join(dir.out.processing, "fa2ids.tsv"),
        output:
            temp(os.path.join(dir.out.processing, "fasta", "circular.done")),
        params:
            seed=SEED,
            circular=CIRCULAR,
            fadir=os.path.abspath(os.path.join(dir.out.processing, "fasta")),
            rotates=ROTATES,
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        script:
            os.path.join(dir.scripts, "circular.py")


fa_stats_input = os.path.join(dir.out.processing, "fasta", "renamed.done")
if PRIMERSEARCH:
    fa_stats_input = os.path.join(dir.out.processing, "fasta", "amplicons.done")

if CIRCULAR:
    fa_stats_input = os.path.join(dir.out.processing, "fasta", "circular.done")


rule get_fasta_stats:
    input:
        fa_stats_input,
    output:
        os.path.join(dir.out.processing, "fasta_stats.tsv"),
    params:
        os.path.join(dir.out.processing, "fasta", "*.fasta")
        if not PRIMERSEARCH
        else os.path.join(dir.out.processing, "primersearch", "*.fasta"),
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
