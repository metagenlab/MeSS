art_args = ""
if CUSTOM_ERR == None:
    art_args += f"-ss {ERROR} "
if CUSTOM_ERR:
    custom_profile = [
        os.path.abspath(file) for file in sorted(glob.glob(f"{CUSTOM_ERR}*"))
    ]
    if PAIRED:
        art_args += f"-1 {custom_profile[0]} -2 {custom_profile[1]} "
    else:
        art_args += f"-1 {custom_profile[0]} "


if PAIRED:
    art_args += f"-p -m {FRAG_LEN} -s {FRAG_SD} "


if BAM:
    art_args += "-sam -M"


fq_prefix = os.path.join(dir.out.short, "{sample}", "{fasta}", "{contig}")
if CIRCULAR:
    fq_prefix = os.path.join(dir.out.short, "{sample}", "{fasta}", "{contig}_{n}")

sam_out = temp(fq_prefix + ".txt")
if BAM:
    sam_out = temp(fq_prefix + ".sam")


fastq_out = [
    temp(fq_prefix + "1.fq"),
    temp(fq_prefix + "2.fq"),
]

if not PAIRED:
    fastq_out = temp(fq_prefix + ".fq")

fasta = os.path.join(dir.out.processing, "split", "{fasta}_{contig}.fna")
if CIRCULAR:
    fasta = os.path.join(
        dir.out.processing, "rotate", "{sample}", "{fasta}_{contig}_{n}.fna"
    )


rule art_illumina:
    input:
        fasta,
    output:
        sam=sam_out,
        fastqs=fastq_out,
    params:
        args=art_args,
        read_len=MEAN_LEN,
        cov=lambda wildcards: get_value("cov_sim", wildcards),
        seed=lambda wildcards: int(get_value("seed", wildcards)),
        prefix=fq_prefix,
    log:
        os.path.join(dir.out.logs, "art", "{sample}", "{fasta}", "{contig}.log")
        if not CIRCULAR
        else os.path.join(
            dir.out.logs, "art", "{sample}", "{fasta}", "{contig}_{n}.log"
        ),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.norm.time,
    conda:
        os.path.join(dir.conda, "art.yml")
    container:
        containers.art
    shell:
        """
        art_illumina -i {input} \\
        -rs {params.seed} -l {params.read_len} \\
        -f {params.cov} -na {params.args} \\
        -o {params.prefix} &> {log}
        touch {output.sam}
        """
