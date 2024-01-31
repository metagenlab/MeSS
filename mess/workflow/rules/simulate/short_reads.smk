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


sam_out = []
if BAM and PAIRED:
    sam_out = temp(os.path.join(dir.out.short, "{sample}", "{fasta}.sam"))
elif BAM and not PAIRED:
    sam_out = temp(os.path.join(dir.out.short, "{sample}", "{fasta}1.sam"))
if not BAM:
    sam_out = temp(os.path.join(dir.out.short, "{sample}", "{fasta}.txt"))

fastq_out = [
    temp(os.path.join(dir.out.short, "{sample}", "{fasta}1.fq")),
    temp(os.path.join(dir.out.short, "{sample}", "{fasta}2.fq")),
]
fq_prefix = os.path.join(dir.out.short, "{sample}", "{fasta}")
if not PAIRED:
    fastq_out = fastq_out[0]
    fq_prefix = os.path.join(dir.out.short, "{sample}", "{fasta}1")

fa = []
if not SKIP_FA_PROC:
    fa = os.path.join(dir.out.fasta, "{fasta}.merged")
else:
    fa = fasta_input


rule art_illumina:
    input:
        fa=fa,
        df=os.path.join(dir.out.base, "cov.tsv"),
    output:
        sam=sam_out,
        fastqs=fastq_out,
    conda:
        os.path.join(dir.env, "art.yml")
    params:
        args=art_args,
        read_len=MEAN_LEN,
        cov=lambda wildcards, input: get_value(input.df, wildcards, "cov_sim"),
        seed=lambda wildcards, input: int(get_value(input.df, wildcards, "seed")),
        prefix=fq_prefix,
    benchmark:
        os.path.join(dir.out.bench, "art", "{sample}", "{fasta}.txt")
    log:
        os.path.join(dir.out.logs, "art", "{sample}", "{fasta}.log"),
    shell:
        """
        art_illumina -i {input.fa}  \\
        -rs {params.seed} -l {params.read_len} \\
        -f {params.cov} -na {params.args} \\
        -o {params.prefix} &> {log}
        touch {output.sam}
        """
