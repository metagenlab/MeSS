"""
Snakefile for simulating reads
"""
import attrmap as ap


# config file
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")


config = ap.AttrMap(config)


wildcard_constraints:
    sample="[^/]+",
    fasta="[^/]+",
    contig="[^/]+",


# functions
include: os.path.join("rules", "preflight", "functions.smk")
# directories
include: os.path.join("rules", "preflight", "directories.smk")


# common options
INPUT = os.path.abspath(str(config.args.input))
OUTPUT = config.args.output
LOG = os.path.join(OUTPUT, "mess.log")
THREADS = config.args.threads
# samples and replicates options
REPLICATES = list(range(1, config.args.replicates + 1))
SAMPLES = parse_samples(INPUT, REPLICATES)
SEED = config.args.seed
REP_SD = config.args.rep_sd


# aggregate samples tables and make replicates
include: os.path.join("rules", "preflight", "samples.smk")


# coverage options
SEQ_TECH = config.args.tech
ASM_SUMMARY = config.args.asm_summary
BASES = config.args.bases
PAIRED = config.args.paired
if not SEQ_TECH == "illumina":
    PAIRED = False
if PAIRED:
    FRAG_LEN = config.args.frag_len
    FRAG_SD = config.args.frag_sd
    PAIRS = [1, 2]
else:
    PAIRS = [1]

DIST = config.args.dist
MU = config.args.mu
SIGMA = config.args.sigma
MEAN_LEN = config.args.mean_len


# calculate coverages
include: os.path.join("rules", "processing", "coverages.smk")


# process fasta options
COMPRESSED = config.args.compressed
SKIP_FA_PROC = config.args.skip_fa_proc
if not SKIP_FA_PROC:

    include: os.path.join("rules", "processing", "fastas.smk")


# simulators options
CUSTOM_ERR = config.args.custom_err
ERROR = config.args.error
BAM = config.args.bam
MIN_LEN = config.args.min_len
MAX_LEN = config.args.max_len
SD_LEN = config.args.sd_len
# simulate reads
if SEQ_TECH == "illumina":

    include: os.path.join("rules", "simulate", "short_reads.smk")

else:
    MODEL = ""
    RATIO = ""
    PASSES = config.args.passes
    ACCURACY = config.args.accuracy
    if SEQ_TECH == "pacbio":
        MODEL = "QSHMM-RSII"
        RATIO = "22:45:33"
    elif SEQ_TECH == "nanopore":
        RATIO = "39:24:36"
        if ERROR == "r10.4":
            MODEL = "QSHMM-ONT-HQ"
        elif ERROR == "r10.3":
            MODEL = "QSHMM-ONT"

    include: os.path.join("rules", "simulate", "long_reads.smk")


# process reads
SHUFFLE = dict(zip(SAMPLES, random.sample(range(1, 100000), len(SAMPLES))))


# reads post-processsing options
SKIP_SHUFFLE = config.args.skip_shuffle


include: os.path.join("rules", "processing", "reads.smk")
include: os.path.join("rules", "preflight", "targets_simulate.smk")


# define the rule all
rule simulate:
    input:
        TargetSimreads,
