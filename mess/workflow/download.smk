"""
Snakefile for downloading assemblies
"""
import attrmap as ap
import os
import glob


# Concatenate Snakemake's own log file with the master log file
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    shell("cat " + current_log + " >> " + LOG)


onsuccess:
    copy_log_file()


onerror:
    copy_log_file()


# config file
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")


config = ap.AttrMap(config)


# args
INPUT = os.path.abspath(str(config.args.input))
OUTPUT = config.args.output
LOG = os.path.join(OUTPUT, "mess.log")
THREADS = config.args.threads
CONDA_PREFIX = config.args.conda_prefix
TAXONKIT = config.args.taxonkit
API_KEY = config.args.api_key
LIMIT = config.args.limit
COMPRESSED = config.args.compressed
SOURCE = config.args.source
INCLUDE = config.args.include
TAXON = config.args.taxon
REFERENCE = config.args.reference
ASM_LVL = config.args.assembly_level
ANNOTATED = config.args.annotated
ATYPICAL = config.args.atypical
MAG = config.args.mag
RANK = config.args.rank
NRANK = config.args.nrank


# functions
include: os.path.join("rules", "preflight", "functions.smk")
# directories
include: os.path.join("rules", "preflight", "directories.smk")
# read sample
include: os.path.join("rules", "preflight", "samples.smk")
# targets
include: os.path.join("rules", "preflight", "targets_download.smk")
# assembly_finder
include: os.path.join("rules", "download", "assembly_finder.smk")


# define the rule all
rule download:
    input:
        TargetDownloads,
