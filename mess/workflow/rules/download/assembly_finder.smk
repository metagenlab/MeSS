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


af_args = ""
if TAXON:
    af_args += "--taxon "
    if LIMIT:
        af_args += f"--limit {LIMIT} "

    if RANK and NRANK:
        af_args += f"--rank {RANK} --nrank {NRANK}"

    if ASM_LVL:
        af_args += f"--assembly-level {ASM_LVL} "
    if INCLUDE:
        f"--include {INCLUDE} "
    if SOURCE:
        f"--source {SOURCE} "
    if REFERENCE:
        f"--reference {REFERENCE} "
    if ANNOTATED:
        f"--annotated {ANNOTATED} "
    if ATYPICAL:
        f"--atypical {ATYPICAL} "
    if MAG:
        f"--mag {MAG}"

else:
    af_args += "--accession "

if API_KEY:
    af_args += f"--api-key {API_KEY} "

if COMPRESSED:
    af_args += f"--compressed {COMPRESSED} "


checkpoint download_assemblies:
    input:
        tsv=os.path.join(dir.out.base, "uniq_entries.tsv"),
    output:
        asm=os.path.join(dir.out.base, "assembly_finder/assembly_summary.tsv"),
        seq=os.path.join(dir.out.base, "assembly_finder/sequence_report.tsv"),
        tax=os.path.join(dir.out.base, "assembly_finder/taxonomy.tsv"),
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
        --no-use-conda \\
        -o {params.out} 2> {log}
        """
