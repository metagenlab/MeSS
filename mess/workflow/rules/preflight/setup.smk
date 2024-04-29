rule aggregate_samples:
    output:
        temp(os.path.join(dir.out.base, "samples.tsv")),
    params:
        INPUT,
    script:
        os.path.join(dir.scripts, "samples.py")


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
        os.path.join(dir.env, "utils.yml")
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
    conda:
        os.path.join(dir.env, "utils.yml")
    shell:
        """
        tar -xzvf {input} -C {params.dir} &> {log}
        """
