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
    art_args += "-sam"


sam_out = temp(os.path.join(dir.out.short, "{sample}", "{fasta}", "{contig}.txt"))
if BAM:
    sam_out = temp(os.path.join(dir.out.short, "{sample}", "{fasta}", "{contig}.sam"))


fastq_out = [
    temp(os.path.join(dir.out.short, "{sample}", "{fasta}", "{contig}1.fq")),
    temp(os.path.join(dir.out.short, "{sample}", "{fasta}", "{contig}2.fq")),
]

fq_prefix = os.path.join(dir.out.short, "{sample}", "{fasta}", "{contig}")

if not PAIRED:
    fastq_out = temp(os.path.join(dir.out.short, "{sample}", "{fasta}", "{contig}.fq"))


rule art_illumina:
    input:
        fa=os.path.join(dir.out.processing, "split", "{fasta}_{contig}.fna"),
        df=get_cov_table,
    output:
        sam=sam_out,
        fastqs=fastq_out,
    params:
        args=art_args,
        read_len=MEAN_LEN,
        cov=lambda wildcards, input: get_value(input.df, wildcards, "cov_sim"),
        seed=lambda wildcards, input: int(get_value(input.df, wildcards, "seed")),
        prefix=fq_prefix,
    benchmark:
        os.path.join(dir.out.bench, "art", "{sample}", "{fasta}", "{contig}.txt")
    log:
        os.path.join(dir.out.logs, "art", "{sample}", "{fasta}", "{contig}.log"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.norm.time,
    conda:
        os.path.join(dir.conda, "art.yml")
    shell:
        """
        art_illumina -i {input.fa} \\
        -rs {params.seed} -l {params.read_len} \\
        -f {params.cov} -na {params.args} \\
        -o {params.prefix} &> {log}
        touch {output.sam}
        """
