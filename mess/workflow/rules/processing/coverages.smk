rule replicates_table:
    input:
        os.path.join(dir.out.base, "samples.tsv"),
    output:
        temp(os.path.join(dir.out.base, "replicates.tsv")),
    params:
        seed=SEED,
        replicates=REPLICATES,
        rep_sd=REP_SD,
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    log:
        os.path.join(dir.out.logs, "tables", "replicates.tsv"),
    script:
        os.path.join(dir.scripts, "replicates.py")


checkpoint calculate_coverage:
    input:
        df=os.path.join(dir.out.base, "replicates.tsv"),
        asm=get_assembly_summary_path,
    output:
        temp(
            os.path.join(dir.out.base, "cov.tsv"),
        ),
    params:
        dist=DIST,
        mu=MU,
        sigma=SIGMA,
        bases=BASES,
        read_len=MEAN_LEN,
        pairing=PAIRED,
        seed=SEED,
    log:
        os.path.join(dir.out.logs, "tables", "cov.tsv"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    script:
        os.path.join(dir.scripts, "calculate_cov.py")
