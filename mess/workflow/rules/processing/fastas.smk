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


checkpoint split_contigs:
    input:
        fa=list_fastas,
        cov=os.path.join(dir.out.processing, "cov.tsv"),
    output:
        tsv=temp(os.path.join(dir.out.processing, "split.tsv")),
        dir=temp(directory(os.path.join(dir.out.base, "split"))),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    log:
        os.path.join(dir.out.logs, "tables", "split.tsv"),
    script:
        os.path.join(dir.scripts, "split_contigs.py")
