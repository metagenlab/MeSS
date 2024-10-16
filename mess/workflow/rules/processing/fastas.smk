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
            fa=os.path.join(dir.out.processing, "split", "{fasta}_{contig}.fna"),
            df=get_cov_table,
        output:
            os.path.join(dir.out.processing, "rotate", "{fasta}_{contig}_{n}.fna"),
        params:
            random_start=get_random_start,
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        run:
            for record in SeqIO.parse(input[0], "fasta"):
                header = record.id + "_" + wildcards.n
                record.name = header
                record.id = header
                record.description = header
                record.seq = (
                    record.seq[params.random_start :]
                    + record.seq[: params.random_start]
                )
                SeqIO.write(record, output[0], "fasta")

