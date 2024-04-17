rule aggregate_samples:
    output:
        temp(os.path.join(dir.out.base, "samples.tsv")),
    params:
        INPUT,
    script:
        os.path.join(dir.scripts, "samples.py")
