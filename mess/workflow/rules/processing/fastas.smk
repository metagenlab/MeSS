if GZIP:

    rule unzip_fasta:
        input:
            fasta_input,
        output:
            temp(os.path.join(dir.out.fasta, "{fasta}.fna")),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        shell:
            "zcat {input} > {output}"


rule rename_headers:
    input:
        os.path.join(dir.out.fasta, "{fasta}.fna") if GZIP else fasta_input,
    output:
        temp(os.path.join(dir.out.fasta, "{fasta}.fasta")),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    shell:
        """
        cut -f 1 -d" " {input} > {output}
        """


checkpoint split_contigs:
    input:
        fa=list_fastas,
        cov=os.path.join(dir.out.base, "cov.tsv"),
    output:
        tsv=temp(os.path.join(dir.out.base, "split.tsv")),
        dir=temp(directory(os.path.join(dir.out.base, "split"))),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    conda:
        os.path.join(dir.env, "fasta.yml")
    log:
        os.path.join(dir.out.logs, "tables", "split.tsv"),
    script:
        os.path.join(dir.scripts, "split_contigs.py")
