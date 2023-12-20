checkpoint calculate_coverage:
    input:
        get_assembly_summary_path,
    output:
        temp(
            os.path.join(dir.out.base, "cov.tsv"),
        ),
    params:
        dist=DIST,
        mu=MU,
        sigma=SIGMA,
        total_bases=TOTAL_BASES,
        read_len=MEAN_LEN,
        pairing=PAIRED,
        rep_sd=REP_SD,
        seed=SEED,
        chunks=CHUNKS,
    log:
        os.path.join(dir.out.logs, "tables", "cov.tsv"),
    benchmark:
        os.path.join(dir.out.bench, "tables", "cov.txt")
    script:
        os.path.join(dir.scripts, "calculate_cov.py")
