import os
import click
from snaketool_utils.cli_utils import echo_click


def snake_base(rel_path):
    """Get the filepath to a Snaketool system file (relative to __main__.py)"""
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    """Read and print the version from the version file"""
    with open(snake_base("mess.VERSION"), "r") as f:
        version = f.readline()
    return version


def print_citation():
    """Read and print the Citation information from the citation file"""
    with open(snake_base("mess.CITATION"), "r") as f:
        for line in f:
            echo_click(line)


def default_to_output(ctx, param, value):
    """Callback for click options; places value in output directory unless specified"""
    if param.default == value:
        return os.path.join(ctx.params["output"], value)
    return value


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the common_options decorator.
    """
    options = [
        click.option(
            "-o",
            "--output",
            help="Output directory",
            type=click.Path(dir_okay=True, writable=True, readable=True),
            default="mess_out",
            show_default=True,
        ),
        click.option(
            "--configfile",
            default="config.yaml",
            show_default=False,
            callback=default_to_output,
            help="Custom config file [default: (outputDir)/config.yaml]",
        ),
        click.option(
            "--threads", help="Number of threads to use", default=1, show_default=True
        ),
        click.option(
            "--profile",
            default=None,
            help="Snakemake profile to use",
            show_default=False,
        ),
        click.option(
            "--use-conda/--no-use-conda",
            default=True,
            help="Use conda for Snakemake rules",
            show_default=True,
        ),
        click.option(
            "--conda-prefix",
            default=snake_base(os.path.join("workflow", "conda")),
            help="Custom conda env directory",
            type=click.Path(),
            show_default=False,
        ),
        click.option(
            "--snake-default",
            multiple=True,
            default=[
                "--printshellcmds",
                "--nolock",
                "--show-failed-logs",
            ],
            help="Customise Snakemake runtime args",
            show_default=True,
        ),
        click.option(
            "--log",
            default="mess.log",
            callback=default_to_output,
            hidden=True,
        ),
        click.option(
            "--system-config",
            default=snake_base(os.path.join("config", "config.yaml")),
            hidden=True,
        ),
        click.argument("snake_args", nargs=-1),
    ]
    for option in reversed(options):
        func = option(func)
    return func


def sim_options(func):
    """Reads simulation command line args
    include simulate args with the sim_options decorator.
    """
    options = [
        click.option(
            "--tech",
            help="sequencing technology",
            type=click.Choice(["illumina", "pacbio", "nanopore"], case_sensitive=False),
            default="illumina",
            show_default=True,
        ),
        click.option(
            "--paired",
            help="illumina reads pairing",
            type=bool,
            default=True,
            show_default=True,
        ),
        click.option(
            "--tool",
            help="choose simulator",
            type=click.Choice(
                ["art", "pbsim3"],
                case_sensitive=False,
            ),
            default="art",
            show_default=True,
        ),
        click.option(
            "--error",
            help="choose simulator error profile",
            type=str,
            default="MSv3",
            show_default=True,
        ),
        click.option(
            "--replicates",
            help="number of replicates per sample",
            type=int,
            default=1,
            show_default=True,
        ),
        click.option(
            "--rep-sd",
            help="standard deviation between replicate coverages (0: no deviation)",
            type=int,
            default=0,
            show_default=True,
        ),
        click.option(
            "--dist",
            help="sample genome coverage values from a distribution",
            type=click.Choice(["none", "lognormal", "even"]),
            default="none",
            show_default=True,
        ),
        click.option(
            "--mu",
            help="mu of lognormal dist",
            type=int,
            default=1,
        ),
        click.option(
            "--sigma",
            help="sigma of lognormal dist",
            type=int,
            default=0,
        ),
        click.option(
            "--bases",
            help="per sample base counts (ignored when coverage is set)",
            type=str,
            default="1G",
            show_default=True,
        ),
        click.option(
            "--accuracy",
            help="mean accuracy for long read sequencing",
            type=float,
            default=0.9,
            show_default=True,
        ),
        click.option(
            "--passes",
            help="number of passes for pacbio multipass sequencing (2 or more for hifi)",
            type=int,
            default=1,
            show_default=True,
        ),
        click.option(
            "--min-len",
            help="minimum read length for long read sequencing",
            type=int,
            default=100,
            show_default=True,
        ),
        click.option(
            "--max-len",
            help="maximum read length for long read sequencing",
            type=int,
            default=300,
            show_default=True,
        ),
        click.option(
            "--sd-len",
            help="standard read length deviation for long read sequencing",
            type=int,
            default=50,
            show_default=True,
        ),
        click.option(
            "--mean-len",
            help="mean read length for long and short read sequencing",
            type=int,
            default=150,
            show_default=True,
        ),
        click.option(
            "--frag-len",
            help="Fragment length for paird short read sequencing",
            type=int,
            default=170,
            show_default=True,
        ),
        click.option(
            "--frag-sd",
            help="Fragment length standard deviation for paird short read sequencing",
            type=int,
            default=17,
            show_default=True,
        ),
        click.option(
            "--bam",
            help="generate gold standard alignment bam files",
            type=bool,
            default=False,
            show_default=True,
        ),
        click.option(
            "--seed",
            help="seed for pseudo random number generator",
            type=int,
            default=42,
            show_default=True,
        ),
    ]

    for option in reversed(options):
        func = option(func)
    return func


def download_options(func):
    """Assembly finder command line args
    include download args with the download_options decorator.
    """
    options = [
        click.option(
            "--db",
            type=click.Choice(["refseq", "genbank"], case_sensitive=False),
            help="download from refseq or genbank",
            default="refseq",
            show_default=True,
        ),
        click.option(
            "--uid",
            help="are inputs UIDs or assembly names",
            type=str,
            default=False,
            show_default=True,
        ),
        click.option(
            "--rc",
            type=str,
            help="select reference, representative or all",
            default="all",
            show_default=True,
        ),
        click.option(
            "--al",
            type=str,
            help="select complete, chromosome, scaffold, contig or all",
            default="complete",
            show_default=True,
        ),
        click.option(
            "--an",
            type=click.Choice(["False", "True"]),
            help="select assemblies with annotation",
            default="False",
            show_default=True,
        ),
        click.option(
            "--rank",
            help="taxonomic rank to filter by assemblies ",
            default="none",
            type=click.Choice(
                [
                    "superkingdom",
                    "phylum",
                    "class",
                    "order",
                    "family",
                    "genus",
                    "species",
                    "none",
                ],
                case_sensitive=False,
            ),
            show_default=True,
        ),
        click.option(
            "--nr",
            help="max number of genome by target rank (example: 1 per species)",
            type=str,
            default="none",
            show_default=True,
        ),
        click.option(
            "--ete_db",
            type=str,
            default=os.path.join(os.environ["HOME"], ".etetoolkit"),
            help="path where to save/find ete taxa.sqlite file",
            show_default=True,
        ),
        click.option(
            "--ncbi_key", help="NCBI key", type=str, required=False, default="none"
        ),
        click.option(
            "--ncbi_email", help="NCBI email", type=str, required=False, default="none"
        ),
        click.option(
            "--nb",
            help="number of assemblies per entry",
            type=str,
            default="all",
            show_default=True,
        ),
        click.option(
            "--suffixes",
            type=str,
            help="suffix of files to download from NCBI's ftp",
            default="assembly_report.txt,genomic.fna.gz",
            show_default=True,
        ),
        click.option(
            "--exclude",
            type=str,
            help="filter to exclude assemblies (example: exclude from metagenomes)",
            default="metagenome",
            show_default=True,
        ),
    ]

    for option in reversed(options):
        func = option(func)
    return func
