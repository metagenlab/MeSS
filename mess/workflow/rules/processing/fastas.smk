
rule unzip_fasta:
    input:
        get_fasta,
    output:
        temp(os.path.join(dir.out.fasta, "{fasta}.fa")),
    benchmark:
        os.path.join(dir.out.bench, "fasta", "unzip", "{fasta}.txt")
    shell:
        "zcat {input} > {output}"


rule rename_headers:
    input:
        os.path.join(dir.out.fasta, "{fasta}.fa"),
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
