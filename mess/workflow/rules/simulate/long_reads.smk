prefix = os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}")
if PASSES > 1:
    pbsim3_out = temp(os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}.sam"))
    rename = f"mv {prefix}_0001.sam {prefix}.sam"
else:
    pbsim3_out = temp(os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}.fq"))
    rename = f"mv {prefix}_0001.fastq {prefix}.fq"


rule pbsim3:
    input:
        fa=os.path.join(dir.out.processing, "split", "{fasta}_{contig}.fna"),
        df=get_cov_table,
    output:
        pbsim3_out,
        temp(os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}.maf")),
        temp(os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}.ref")),
    params:
        model=os.path.join(MODEL_PATH, f"{MODEL}.model"),
        ratio=RATIO,
        meanlen=MEAN_LEN,
        lensd=SD_LEN,
        minlen=MIN_LEN,
        maxlen=MAX_LEN,
        passes=PASSES,
        accuracy=ACCURACY,
        cov=lambda wildcards, input: get_value(input.df, wildcards, "cov_sim"),
        seed=lambda wildcards, input: int(get_value(input.df, wildcards, "seed")),
        prefix=prefix,
        reads_rename=rename,
    benchmark:
        os.path.join(dir.out.bench, "pbsim3", "{sample}", "{fasta}", "{contig}.txt")
    log:
        os.path.join(dir.out.logs, "pbsim3", "{sample}", "{fasta}", "{contig}.log"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.norm.time,
    conda:
        os.path.join(dir.conda, "pbsim3.yml")
    shell:
        """
        pbsim --strategy wgs --method qshmm \\
        --id-prefix {wildcards.contig} \\
        --difference-ratio {params.ratio} \\
        --length-min {params.minlen} \\
        --length-max {params.maxlen} \\
        --length-mean {params.meanlen} \\
        --length-sd {params.lensd} \\
        --prefix {params.prefix} \\
        --seed {params.seed} \\
        --qshmm {params.model} \\
        --pass-num {params.passes} \\
        --accuracy-mean {params.accuracy} \\
        --depth {params.cov} --genome {input.fa} &> {log}
        
        mv {params.prefix}_0001.maf {params.prefix}.maf
        mv {params.prefix}_0001.ref {params.prefix}.ref
        {params.reads_rename}
        """
