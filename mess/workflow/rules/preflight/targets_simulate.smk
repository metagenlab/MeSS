"""
All simulated reads are declared here
"""


def list_reads(wildcards):
    if SEQ_TECH == "illumina":
        reads = expand(
            os.path.join(dir.out.fastq, "{sample}_R{p}.fq.gz"),
            sample=SAMPLES,
            p=PAIRS,
        )
    else:
        reads = expand(
            os.path.join(dir.out.fastq, "{sample}.fq.gz"),
            sample=SAMPLES,
        )

    if BAM:
        bams = expand(
            os.path.join(dir.out.bam, "{sample}.{bam}"),
            sample=SAMPLES,
            bam=["bam", "bam.bai"],
        )
        reads.append(bams)
    return reads


TargetSimreads = [list_reads, os.path.join(dir.out.base, "cov.tsv")]
