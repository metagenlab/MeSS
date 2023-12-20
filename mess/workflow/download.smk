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


# variables
INPUT = os.path.abspath(str(config.args.input))
OUTPUT = config.args.output
LOG = os.path.join(OUTPUT, "mess.log")
THREADS = config.args.threads
NCBI_EMAIL = config.args.ncbi_email
NCBI_KEY = config.args.ncbi_key


# functions
include: os.path.join("rules", "preflight", "functions.smk")
# directories
include: os.path.join("rules", "preflight", "directories.smk")
# read sample
include: os.path.join("rules", "preflight", "read_samples.smk")
# targets
include: os.path.join("rules", "preflight", "targets_download.smk")
# assembly_finder
include: os.path.join("rules", "download", "assembly_finder.smk")


# define the rule all
rule download:
    input:
        TargetDownloads,
