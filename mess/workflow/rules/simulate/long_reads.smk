prefix = os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}")
fasta = os.path.join(dir.out.processing, "split", "{fasta}_{contig}.fna")
if CIRCULAR:
    prefix = os.path.join(dir.out.long, "{sample}", "{fasta}", "{contig}_{n}")
    fasta = os.path.join(
        dir.out.processing, "rotate", "{sample}", "{fasta}_{contig}_{n}.fna"
    )
id_prefix = os.path.basename(prefix)


if PASSES > 1:
    pbsim3_out = temp(prefix + ".sam")
    rename = f"mv {prefix}_0001.sam {prefix}.sam"
else:
    pbsim3_out = temp(prefix + ".fq")
    rename = f"mv {prefix}_0001.fastq {prefix}.fq"


rule pbsim3:
    input:
        df=get_cov_df,
        fa=fasta,
    output:
        pbsim3_out,
        temp(prefix + ".maf"),
        temp(prefix + ".ref"),
    params:
        model=os.path.join(MODEL_PATH, f"{MODEL}.model"),
        ratio=RATIO,
        meanlen=MEAN_LEN,
        lensd=SD_LEN,
        minlen=MIN_LEN,
        maxlen=MAX_LEN,
        passes=PASSES,
        accuracy=ACCURACY,
        cov=lambda wildcards: get_value("cov_sim", wildcards),
        seed=lambda wildcards: int(get_value("seed", wildcards)),
        prefix=prefix,
        id_prefix=id_prefix,
        rename=rename,
    log:
        os.path.join(dir.out.logs, "pbsim3", "{sample}", "{fasta}", "{contig}.log")
        if not CIRCULAR
        else os.path.join(
            dir.out.logs, "pbsim3", "{sample}", "{fasta}", "{contig}_{n}.log"
        ),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.norm.time,
    conda:
        os.path.join(dir.conda, "pbsim3.yml")
    container:
        containers.pbsim3
    shell:
        """
        pbsim --strategy wgs --method qshmm \\
        --id-prefix {params.id_prefix} \\
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
        --depth {params.cov} \\
        --genome {input.fa} &> {log}
        
        mv {params.prefix}_0001.maf {params.prefix}.maf
        mv {params.prefix}_0001.ref {params.prefix}.ref
        {params.rename}
        """
