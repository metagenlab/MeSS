"""
Snakefile for downloading assemblies
"""


include: os.path.join("rules", "preflight", "config.smk")


# assembly_finder options
REFERENCE = config.args.reference
LIMIT = config.args.limit
TAXON = config.args.taxon
AF_ARGS = config.args.af_args


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
