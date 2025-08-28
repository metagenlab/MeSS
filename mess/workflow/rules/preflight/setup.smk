rule aggregate_samples:
    output:
        temp(os.path.join(dir.out.base, "samples.tsv")),
    params:
        INPUT,
    script:
        os.path.join(dir.scripts, "samples.py")


if not os.path.exists(os.path.join(TAXONKIT, "names.dmp")):

    rule download_taxdump:
        output:
            os.path.join(TAXONKIT, "taxdump.tar.gz"),
        log:
            os.path.join(dir.out.logs, "curl.log"),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        conda:
            os.path.join(dir.conda, "curl.yml")
        container:
            containers.curl
        shell:
            """
            curl https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz \\
            -o {output} 2> {log}
            """

    rule decompress_taxdump:
        input:
            os.path.join(TAXONKIT, "taxdump.tar.gz"),
        output:
            os.path.join(TAXONKIT, "names.dmp"),
            os.path.join(TAXONKIT, "nodes.dmp"),
            os.path.join(TAXONKIT, "delnodes.dmp"),
            os.path.join(TAXONKIT, "merged.dmp"),
        params:
            dir=TAXONKIT,
        log:
            os.path.join(dir.out.logs, "taxdump.log"),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        shell:
            """
            tar -xzvf {input} -C {params.dir} &> {log}
            """


RANKS = config.args.ranks
if isinstance(config.args.custom_tax, str):
    RANKS = ",".join(list(custom_tax_df.columns[1:]))
    TAXONKIT = os.path.join(TAXONKIT, "custom")

    rule create_custom_taxdump:
        input:
            os.path.abspath(config.args.custom_tax),
        output:
            os.path.join(TAXONKIT, "names.dmp"),
            os.path.join(TAXONKIT, "nodes.dmp"),
            os.path.join(TAXONKIT, "delnodes.dmp"),
            os.path.join(TAXONKIT, "merged.dmp"),
            os.path.join(TAXONKIT, "taxid.map"),
        params:
            outdir=TAXONKIT,
            ranks=RANKS,
        log:
            os.path.join(dir.out.logs, "custom-taxdump.log"),
        resources:
            mem_mb=config.resources.sml.mem,
            mem=str(config.resources.sml.mem) + "MB",
            time=config.resources.sml.time,
        conda:
            os.path.join(dir.conda, "taxonkit.yml")
        container:
            containers.taxonkit
        shell:
            """
            taxonkit create-taxdump -A 1 \\
            -R {params.ranks} \\
            {input} \\
            -O {params.outdir} \\
            2> {log}
            """


rule get_primers_table:
    output:
        os.path.join(dir.out.processing, "primers.txt"),
    params:
        fw=config.args.fw,
        rv=config.args.rv,
    run:
        fw = params.fw.split(",")
        rv = params.rv.split(",")
        names = [f"primer_{n+1}" for n in range(len(fw))]
        pd.DataFrame(list(zip(names, fw, rv))).to_csv(
            output[0], sep="\t", index=False, header=False
        )
