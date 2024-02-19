rule get_unique_entries:
    input:
        os.path.join(dir.out.base, "samples.tsv"),
    output:
        temp(os.path.join(dir.out.base, "uniq_entries.tsv")),
    params:
        nb=NB,
    run:
        df = pd.read_csv(input[0], sep="\t")
        try:
            df[["entry", "nb"]].drop_duplicates().to_csv(
                output[0], sep="\t", index=None
            )
        except KeyError:
            df["nb"] = [params.nb] * len(df)
            df[["entry", "nb"]].drop_duplicates().to_csv(
                output[0], sep="\t", index=None
            )


checkpoint download_assemblies:
    input:
        os.path.join(dir.out.base, "uniq_entries.tsv"),
    output:
        asm=os.path.join(dir.out.base, "download/assembly_summary.tsv"),
        version=os.path.join(dir.out.versions, "assembly_finder.version"),
    params:
        ncbi_key=NCBI_KEY,
        ncbi_email=NCBI_EMAIL,
        db=DB,
        uid=UID,
        alvl=ASM_LVL,
        refc=REF_CAT,
        annot=ANNOT,
        excl=EXCLUDE,
        rank=RANK,
        nr=NB_RANK,
        ete=ETE_DB,
        out=os.path.join(dir.out.base, "download"),
    benchmark:
        os.path.join(dir.out.bench, "download", "download.txt")
    log:
        os.path.join(dir.out.logs, "download", "download.log"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.med.cpu
    conda:
        os.path.join(dir.env, "assembly_finder.yml")
    shell:
        """
        assembly_finder \\
        -i {input} \\
        -t {threads} \\
        -nk {params.ncbi_key} \\
        -ne {params.ncbi_email} \\
        -db {params.db} \\
        -id {params.uid} \\
        -al {params.alvl} \\
        -rc {params.refc} \\
        -an {params.annot} \\
        -ex {params.excl} \\
        -r {params.rank} \\
        -nr {params.nr} \\
        -et {params.ete} \\
        -o {params.out} \\
        --nolock 2> {log}
        assembly_finder -v > {output.version}
        """
