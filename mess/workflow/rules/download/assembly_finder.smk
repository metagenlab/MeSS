rule get_unique_entries:
    input:
        os.path.join(dir.out.base, "samples.tsv"),
    output:
        temp(os.path.join(dir.out.base, "uniq_entries.tsv")),
    params:
        nb=LIMIT,
        taxon=TAXON,
    run:
        df = pd.read_csv(input[0], sep="\t")
        if params.taxon:
            if "nb" not in df.columns:
                df["nb"] = [params.nb] * len(df)
            df[["taxon", "nb"]].drop_duplicates().to_csv(
                output[0], sep="\t", index=False
            )
        else:
            df[["accession"]].drop_duplicates().to_csv(
                output[0], sep="\t", index=False, header=False
            )


af_args = AF_ARGS + " "
if REFERENCE:
    af_args += "--reference "
if TAXON:
    af_args += "--taxon "
    if LIMIT:
        af_args += f"--limit {LIMIT} "

else:
    af_args += "--accession "

af_args = af_args.strip()


checkpoint download_assemblies:
    input:
        tsv=os.path.join(dir.out.base, "uniq_entries.tsv"),
    output:
        asm=os.path.join(dir.out.base, "assembly_finder", "assembly_summary.tsv"),
        tax=os.path.join(dir.out.base, "assembly_finder", "taxonomy.tsv"),
    params:
        args=af_args,
        taxonkit=TAXONKIT,
        out=os.path.join(dir.out.base, "assembly_finder"),
    log:
        os.path.join(dir.out.logs, "assembly_finder.log"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.norm.cpu
    conda:
        os.path.join(dir.conda, "assembly_finder.yml")
    container:
        containers.assembly_finder
    shell:
        """
        assembly_finder -i {input.tsv} \\
        --taxonkit {params.taxonkit} \\
        --threads {threads} \\
        {params.args} \\
        -o {params.out} 2> {log}
        """
