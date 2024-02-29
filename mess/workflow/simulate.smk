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
    p="[1-2]+",


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
PASSES = config.args.passes
ACCURACY = config.args.accuracy
MODEL = config.args.model
RATIO = config.args.ratio

if SEQ_TECH == "illumina":

    include: os.path.join("rules", "simulate", "short_reads.smk")

else:

    include: os.path.join("rules", "simulate", "long_reads.smk")


# reads post-processsing options
SHUFFLE = dict(zip(SAMPLES, random.sample(range(1, 100000), len(SAMPLES))))
SKIP_SHUFFLE = config.args.skip_shuffle


include: os.path.join("rules", "processing", "reads.smk")
include: os.path.join("rules", "preflight", "targets_simulate.smk")


# define the rule all
rule simulate:
    input:
        TargetSimreads,
