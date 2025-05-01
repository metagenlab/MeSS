"""
Snakefile for downloading assemblies
"""


include: os.path.join("rules", "preflight", "config.smk")


# assembly_finder options
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


# assembly_finder rules
include: os.path.join("rules", "download", "assembly_finder.smk")


# assembly_finder outputs
TargetDownloads = [
    os.path.join(dir.out.base, "uniq_entries.tsv"),
    os.path.join(dir.out.base, "assembly_finder", "assembly_summary.tsv"),
    os.path.join(dir.out.base, "assembly_finder", "taxonomy.tsv"),
]


# define the rule all
rule download:
    input:
        TargetDownloads,
