"""
Snakefile for downloading assemblies
"""
import attrmap as ap


# config file
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")


config = ap.AttrMap(config)

# variables


# functions
include: os.path.join("rules", "preflight", "functions.smk")
# directories
include: os.path.join("rules", "preflight", "directories.smk")


# variables
INPUT = os.path.abspath(str(config.args.input))
NCBI_EMAIL = config.args.ncbi_email
NCBI_KEY = config.args.ncbi_key


# read sample
include: os.path.join("rules", "preflight", "read_samples.smk")
# targets
include: os.path.join("rules", "preflight", "target_downloads.smk")
# assembly_finder
include: os.path.join("rules", "download", "assembly_finder.smk")


# define the rule all
rule download:
    input:
        TargetDownloads,
