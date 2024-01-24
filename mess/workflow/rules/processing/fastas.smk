if COMPRESSED:

    rule unzip_fasta:
        input:
            fasta_input,
        output:
            temp(os.path.join(dir.out.fasta, "{fasta}.fa")),
        benchmark:
            os.path.join(dir.out.bench, "fasta", "unzip", "{fasta}.txt")
        shell:
            "zcat {input} > {output}"


rule rename_headers:
    input:
        os.path.join(dir.out.fasta, "{fasta}.fa") if COMPRESSED else fasta_input,
    output:
        temp(os.path.join(dir.out.fasta, "{fasta}.renamed")),
    benchmark:
        os.path.join(dir.out.bench, "fasta", "rename", "{fasta}.txt")
    shell:
        """
        cut -f 1 -d" " {input} > {output}
        """


rule merge_contigs:
    input:
        os.path.join(dir.out.fasta, "{fasta}.renamed"),
    output:
        temp(os.path.join(dir.out.fasta, "{fasta}.merged")),
    benchmark:
        os.path.join(dir.out.bench, "fasta", "merge", "{fasta}.txt")
    script:
        os.path.join(dir.scripts, "merge_contigs.py")


checkpoint split_contigs:
    input:
        lambda wildcards: list_fastas(wildcards, contigs=True),
    output:
        temp(os.path.join(dir.out.fasta, "split", "split.tsv")),
    params:
        os.path.join(dir.out.fasta, "split"),
    script:
        os.path.join(dir.scripts, "split_contigs.py")


rule remove_split_fasta:
    input:
        os.path.join(dir.out.fasta, "split", "split.tsv"),
    output:
        temp(os.path.join(dir.out.base, "clean.txt")),
    params:
        os.path.join(dir.out.fasta, "split"),
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
