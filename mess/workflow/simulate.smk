"""
Snakefile for simulating reads
"""
import attrmap as ap


# config file
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")


config = ap.AttrMap(config)


# functions
include: os.path.join("rules", "preflight", "functions.smk")
# directories
include: os.path.join("rules", "preflight", "directories.smk")


# common options
INPUT = os.path.abspath(str(config.args.input))
OUTPUT = config.args.output
LOG = os.path.join(OUTPUT, "mess.log")
THREADS = config.args.threads


# read sample table
include: os.path.join("rules", "preflight", "read_samples.smk")


# replicates options
REPLICATES = list(range(1, config.args.replicates + 1))
REP_SD = config.args.rep_sd
CHUNKS = [f"{x+1:03}" for x in range(config.args.chunks)]


# make replicates table
include: os.path.join("rules", "preflight", "get_replicates.smk")
# process fasta
include: os.path.join("rules", "processing", "fastas.smk")


# coverage options
TOTAL_BASES = config.args.total_bases
PAIRED = config.args.paired
if PAIRED:
    PAIRS = [1, 2]
else:
    PAIRS = [1]

DIST = config.args.dist
MU = config.args.mu
SIGMA = config.args.sigma
MEAN_LEN = config.args.mean_len


# calculate coverages
include: os.path.join("rules", "processing", "coverages.smk")


# simulators options
SEQ_TECH = config.args.seq_tech
BAM = config.args.bam
PASSES = config.args.passes
ACCURACY = config.args.accuracy
MIN_LEN = config.args.min_len
MAX_LEN = config.args.max_len
SD_LEN = config.args.sd_len
# simulate reads
if SEQ_TECH == "illumina":

    include: os.path.join("rules", "simulate", "short_reads.smk")

else:

    include: os.path.join("rules", "simulate", "long_reads.smk")


# process reads
SHUFFLE = dict(zip(SAMPLES, random.sample(range(1, 100000), len(SAMPLES))))


include: os.path.join("rules", "processing", "reads.smk")


# define the rule all
rule simulate:
    input:
        TargetSimreads,
