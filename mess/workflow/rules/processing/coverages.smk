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
    script:
        os.path.join(dir.scripts, "replicates.py")


checkpoint calculate_genome_coverages:
    input:
        rep=os.path.join(dir.out.base, "replicates.tsv"),
        asm=get_asm_summary,
    output:
        os.path.join(dir.out.base, "coverages.tsv"),
    params:
        fa_dir=FASTA_DIR,
        fa_path=FASTA_PATH,
        dist=DIST,
        mu=MU,
        sigma=SIGMA,
        bases=BASES,
        read_len=MEAN_LEN,
        pairing=PAIRED,
        seed=SEED,
        amplicons=PRIMERSEARCH,
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    script:
        os.path.join(dir.scripts, "calculate_cov.py")
