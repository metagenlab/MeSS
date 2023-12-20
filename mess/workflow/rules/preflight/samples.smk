rule samples_table:
    output:
        temp(os.path.join(dir.out.base, "samples.tsv")),
    params:
        input=INPUT,
        seed=SEED,
        chunks=CHUNKS,
        replicates=REPLICATES,
        rep_sd=REP_SD,
    script:
        os.path.join(dir.scripts, "samples.py")
