include: "./functions.smk"


rule samples_table:
    input:
        f"{outdir}/cov.tsv",


rule get_samples_and_replicates:
    input:
        config["input_path"],
    output:
        temp(f"{outdir}/samples.tsv"),
    params:
        rep=config["replicates"],
        rep_sd=config["rep_sd"],
        seed=config["seed"],
    script:
        "get_rep_table.py"


checkpoint calculate_coverage:
    input:
        entry=f"{outdir}/samples.tsv",
        asm=lambda wildcards: get_assembly_summary_path(wildcards),
    output:
        temp(f"{outdir}/cov.tsv"),
    params:
        dist=config["dist"],
        mu=config["mu"],
        sigma=config["sigma"],
        total_bases=config["total_bases"],
        read_len=config["read_len"],
        pairing=config["paired"],
        rep_sd=config["rep_sd"],
        seed=config["seed"],
    log:
        f"{outdir}/logs/tables/cov.tsv",
    script:
        "calculate_cov.py"
