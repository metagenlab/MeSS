"""
Snakefile for simulating reads
"""


include: os.path.join("rules", "preflight", "config.smk")


# samples and replicates options
REPLICATES = list(range(1, config.args.replicates + 1))
SAMPLES = parse_samples(INPUT, REPLICATES)
SEED = config.args.seed
REP_SD = config.args.rep_sd


# fasta paths options
FASTA_DIR = config.args.fasta
if "path" in tsv_df.columns:
    FASTA_PATH = True
    if config.args.sdm == "apptainer":
        fa_dirs = set(os.path.abspath(os.path.dirname(p)) for p in tsv_df["path"])
        for fa_dir in fa_dirs:
            workflow.deployment_settings.apptainer_args += f" -B {fa_dir}:{fa_dir}"

else:
    FASTA_PATH = False

# coverage options
SEQ_TECH = config.args.tech
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


# fasta processing options
ROTATE = config.args.rotate

if ROTATE > 1 or "rotate" in tsv_df.columns:
    CIRCULAR = True
else:
    CIRCULAR = False


include: os.path.join("rules", "processing", "fastas.smk")


# simulators options
CUSTOM_ERR = config.args.custom_err
ERROR = config.args.error
BAM = config.args.bam
ERRFREE = config.args.errfree
MIN_LEN = config.args.min_len
MAX_LEN = config.args.max_len
SD_LEN = config.args.sd_len
PASSES = config.args.passes
ACCURACY = config.args.accuracy
MODEL_PATH = config.args.model_path
MODEL = config.args.model
RATIO = config.args.ratio

# taxonomic profile options
TAX = config.args.tax
CUSTOM_TAX = config.args.custom_tax
if CUSTOM_TAX:
    TAX = True

if SEQ_TECH == "illumina":

    include: os.path.join("rules", "simulate", "short_reads.smk")

else:

    include: os.path.join("rules", "simulate", "long_reads.smk")


# reads post-processsing options
random.seed(SEED)
SHUFFLE = dict(zip(SAMPLES, random.sample(range(1, 100000), len(SAMPLES))))
SKIP_SHUFFLE = config.args.skip_shuffle



include: os.path.join("rules", "processing", "reads.smk")


# simulate outputs
TargetSimreads = [list_reads, os.path.join(dir.out.base, "cleanup.done")]


# define the rule all
rule simulate:
    input:
        TargetSimreads,
