if COMPRESSED:

    rule unzip_fasta:
        input:
            fasta_input,
        output:
            temp(os.path.join(dir.out.fasta, "{fasta}.fa")),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        shell:
            "zcat {input} > {output}"


rule rename_headers:
    input:
        os.path.join(dir.out.fasta, "{fasta}.fa") if COMPRESSED else fasta_input,
    output:
        temp(os.path.join(dir.out.fasta, "{fasta}.renamed")),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    shell:
        """
        cut -f 1 -d" " {input} > {output}
        """


rule merge_contigs:
    input:
        os.path.join(dir.out.fasta, "{fasta}.renamed"),
    output:
        temp(os.path.join(dir.out.fasta, "{fasta}.merged")),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    conda:
        os.path.join(dir.env, "fasta.yml")
    script:
        os.path.join(dir.scripts, "merge_contigs.py")


checkpoint split_contigs:
    input:
        lambda wildcards: list_fastas(wildcards, contigs=True),
    output:
        temp(os.path.join(dir.out.fasta, "split", "split.tsv")),
    params:
        os.path.join(dir.out.fasta, "split"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    conda:
        os.path.join(dir.env, "fasta.yml")
    script:
        os.path.join(dir.scripts, "split_contigs.py")


rule remove_split_fasta:
    input:
        os.path.join(dir.out.fasta, "split", "split.tsv"),
    output:
        temp(os.path.join(dir.out.base, "clean.txt")),
    params:
        os.path.join(dir.out.fasta, "split"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    shell:
        """
        rm -r {params}
        touch {output}
        """


# rule chunk_contigs:
#     input:
#         get_contigs,
#     output:
#         expand(
#             os.path.join(
#                 dir.out.fasta,
#                 "{{fasta}}",
#                 "chunks",
#                 "chunk{chunk}.fa.gz",
#             ),
#             chunk=CHUNKS,
#         ),
#     params:
#         splits=THREADS,
#         out=os.path.join(dir.out.fasta, "{fasta}", "chunks"),
#     shell:
#         """
#         kmcp utils split-genomes \\
#         -m 1 -k 1 --split-number {params.splits} \\
#         --split-overlap 0 {input.fa} -O {params.out}
#         """
# rule decompress_contigs:
#     input:
#         os.path.join(
#             dir.out.fasta,
#             "{fasta}",
#             "chunks",
#             "chunk{chunk}.fa.gz",
#         ),
#     output:
#         os.path.join(
#             dir.out.fasta,
#             "{fasta}",
#             "chunks",
#             "chunk{chunk}.fa",
#         ),
#     shell:
#         "gunzip {input}"
