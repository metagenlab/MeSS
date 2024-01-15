rule get_unique_entries:
    input:
        os.path.join(dir.out.base, "samples.tsv"),
    output:
        temp(os.path.join(dir.out.base, "uniq_entries.tsv")),
    run:
        df = pd.read_csv(input[0], sep="\t")
        df[["entry", "nb"]].drop_duplicates().to_csv(output[0], sep="\t", index=None)


checkpoint download_assemblies:
    input:
        os.path.join(dir.out.base, "uniq_entries.tsv"),
    output:
        asm=os.path.join(dir.out.base, "download/assembly_summary.tsv"),
        version=os.path.join(dir.out.versions, "assembly_finder.version"),
    conda:
        os.path.join(dir.env, "assembly_finder.yml")
    params:
        ncbi_key=NCBI_KEY,
        ncbi_email=NCBI_EMAIL,
        nb=NB,
        db=DB,
        uid=UID,
        alvl=ASM_LVL,
        annot=ANNOT,
        excl=EXCLUDE,
        rank=RANK,
        nr=NB_RANK,
        ete=ETE_DB,
        out=os.path.join(dir.out.base, "download"),
    threads: config.resources.med.cpu
    benchmark:
        os.path.join(dir.out.bench, "download", "download.txt")
    log:
        os.path.join(dir.out.logs, "download", "download.log"),
    shell:
        """
        assembly_finder \
        -i {input} \
        -nk {params.ncbi_key} \
        -ne {params.ncbi_email} \
        -nb {params.nb} \
        -db {params.db} \
        -id {params.uid} \
        -al {params.alvl} \
        -an {params.annot} \
        -ex {params.excl} \
        -r {params.rank} \
        -nr {params.nr} \
        -et {params.ete} \
        -o {params.out} \
        --nolock 2> {log}
        assembly_finder -v > {output.version}
        """
