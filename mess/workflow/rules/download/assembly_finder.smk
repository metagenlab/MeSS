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
                output[0], sep="\t", index=None
            )
        else:
            df[["accession"]].drop_duplicates().to_csv(
                output[0], sep="\t", index=None, header=False
            )


af_args = ""
if TAXON:
    af_args += "--taxon "
else:
    af_args = "--accession "

if LIMIT:
    af_args += f"--limit {LIMIT} "

if RANK and NRANK:
    af_args = f"--rank {RANK} --nrank {NRANK}"

if ASM_LVL:
    af_args = f"--assembly-level {ASM_LVL} "

if API_KEY:
    af_args += f"--api-key {API_KEY} "


checkpoint download_assemblies:
    input:
        os.path.join(dir.out.base, "uniq_entries.tsv"),
    output:
        asm=os.path.join(dir.out.base, "assembly_finder/assembly_summary.tsv"),
        seq=os.path.join(dir.out.base, "assembly_finder/sequence_report.tsv"),
        tax=os.path.join(dir.out.base, "assembly_finder/taxonomy.tsv"),
        version=os.path.join(dir.out.versions, "assembly_finder.version"),
    params:
        args=af_args,
        comp=COMPRESSED,
        source=SOURCE,
        incl=INCLUDE,
        ref=REFERENCE,
        annot=ANNOTATED,
        atyp=ATYPICAL,
        mag=MAG,
        out=os.path.join(dir.out.base, "assembly_finder"),
    benchmark:
        os.path.join(dir.out.bench, "assembly_finder", "download.txt")
    log:
        os.path.join(dir.out.logs, "assembly_finder", "download.log"),
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
        {params.args} \\
        --compressed {params.comp} \\
        --include {params.incl} \\
        --source {params.source} \\
        --reference {params.ref} \\
        --annotated {params.annot} \\
        --atypical {params.atyp} \\
        --mag {params.mag} \\
        -o {params.out} 2> {log}
        assembly_finder -v > {output.version}
        """
