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
        os.path.join(dir.env, "seqkit.yml")
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
            os.path.join(dir.env, "seqkit.yml")
        shell:
            """
            seqkit stats -T -j {threads} {params.path} > {output} 2> {log}
            """


checkpoint split_contigs:
    input:
        fa=list_fastas,
        cov=os.path.join(dir.out.processing, "cov.tsv"),
    output:
        tsv=os.path.join(dir.out.processing, "split.tsv"),
        dir=temp(directory(os.path.join(dir.out.base, "split"))),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    script:
        os.path.join(dir.scripts, "split_contigs.py")
