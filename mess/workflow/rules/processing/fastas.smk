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
        seqkit seq -i {input} > {output}
        """


if FASTA and not ASM_SUMMARY:

    rule get_fasta_stats:
        input:
            FASTA,
        output:
            temp(os.path.join(dir.out.processing, "seqkit_stats.tsv")),
        params:
            path=os.path.join(FASTA, "*"),
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
            seqkit stats -T -j {threads} {params.path} > {output} 2> {log}
            """


checkpoint split_contigs:
    input:
        fa=list_fastas,
        cov=os.path.join(dir.out.base, "coverages.tsv"),
    output:
        tsv=os.path.join(dir.out.processing, "cov.tsv"),
        dir=directory(os.path.join(dir.out.processing, "split")),
    params:
        rotate=ROTATE,
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    script:
        os.path.join(dir.scripts, "split_contigs.py")


if ROTATE > 1:

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
            seqkit restart -i {params} {input} | \\
            seqkit replace -p .+ -r {wildcards.contig}_{wildcards.n} > {output}
            """
