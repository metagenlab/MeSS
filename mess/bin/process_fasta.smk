
rule decompress_assemblies:
    input:
        f"{fasta_dir}/{{fasta}}.fna.gz",
    output:
        temp(f"{outdir}/fasta/{{fasta}}.fa"),
    shell:
        "zcat {input} > {output}"


rule rename_headers:
    input:
        f"{outdir}/fasta/{{fasta}}.fa",
    output:
        temp(f"{outdir}/fasta/{{fasta}}.renamed"),
    shell:
        """
        cut -f 1 -d" " {input} > {output}
        """


rule merge_contigs:
    input:
        f"{outdir}/fasta/{{fasta}}.renamed",
    output:
        temp(f"{outdir}/fasta/{{fasta}}.merged"),
    script:
        "merge_contigs.py"
