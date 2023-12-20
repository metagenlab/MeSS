art_args = ""
if PAIRED:
    art_args += f"-p -m {MAX_LEN} -s {SD_LEN}"
if BAM:
    art_args += " -sam"


sam_out = []
if BAM and PAIRED:
    sam_out = temp(os.path.join(dir.out.fastq, "{sample}", "{fasta}.sam"))
elif BAM and not PAIRED:
    sam_out = temp(os.path.join(dir.out.fastq, "{sample}", "{fasta}1.sam"))
if not BAM:
    sam_out = temp(os.path.join(dir.out.fastq, "{sample}", "{fasta}.txt"))

fastq_out = [
    temp(os.path.join(dir.out.fastq, "{sample}", "{fasta}1.fq")),
    temp(os.path.join(dir.out.fastq, "{sample}", "{fasta}2.fq")),
]
fq_prefix = os.path.join(dir.out.fastq, "{sample}", "{fasta}")
if not PAIRED:
    fastq_out = fastq_out[0]
    fq_prefix = os.path.join(dir.out.fastq, "{sample}", "{fasta}1")


rule art_illumina:
    input:
        fasta=os.path.join(dir.out.fasta, "{fasta}.merged"),
        df=os.path.join(dir.out.base, "cov.tsv"),
    output:
        sam=sam_out,
        fastqs=fastq_out,
    conda:
        os.path.join(dir.env, "art.yml")
    params:
        args=art_args,
        system=PROFILE,
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
        art_illumina -ss {params.system} \\
        -i {input.fasta} -rs {params.seed} \\
        -l {params.read_len} -f {params.cov} \\
        -na {params.args} -o {params.prefix} &> {log}
        touch {output.sam}
        """
