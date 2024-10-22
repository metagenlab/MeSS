import os
import rich_click as click
import yaml
import sys
import subprocess
from snaketool_utils.cli_utils import (
    msg,
    read_config,
    update_config,
    copy_config,
    echo_click,
    msg_box,
)

workflow_basedir = os.path.join(os.path.dirname(os.path.realpath(__file__)))


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


def get_deployment_prefix(ctx, param, value):
    """Callback for click options; joins software deployment name and prefix dir"""
    if not value:
        return os.path.join(os.getcwd(), ".snakemake", ctx.params["sdm"])
    return value


def run_snakemake(
    configfile=None,
    system_config=None,
    snakefile_path=None,
    merge_config=None,
    threads=1,
    sdm=None,
    prefix=None,
    snake_default=None,
    snake_args=[],
    profile=None,
    workflow_profile=None,
    system_workflow_profile=None,
    log=None,
    **kwargs,
):
    """Run a Snakefile!

    Args:
        configfile (str): Filepath of config file to pass with --configfile
        system_config (str): Filepath of system config to copy if configfile not present
        snakefile_path (str): Filepath of Snakefile
        merge_config (dict): Config values to merge with your config file
        threads (int): Number of local threads to request
        deployment (str): Snakemake's --sdm option
        prefix (str): Filepath for Snakemake's --sdm option
        snake_default (list): Snakemake args to pass to Snakemake
        snake_args (list): Additional args to pass to Snakemake
        profile (str): Name of Snakemake profile
        workflow_profile (str): Name of Snakemake workflow-profile
        system_workflow_profile (str): Filepath of system workflow-profile config.yaml to copy if not present
        log (str): Log file for writing STDERR
        **kwargs:

    Returns (int): Exit code
    """

    snake_command = ["snakemake", "-s", snakefile_path]

    # if using a configfile
    if configfile:
        # copy sys default config if not present
        copy_config(configfile, system_config=system_config, log=log)

        if merge_config:
            update_config(in_config=configfile, merge=merge_config, log=log)

        snake_command += ["--configfile", configfile]

        # display the runtime configuration
        snake_config = read_config(configfile)
        msg_box(
            "Runtime config",
            errmsg=yaml.dump(snake_config, Dumper=yaml.Dumper),
            log=log,
        )

    # add threads
    if "--profile" not in snake_args and profile is None:
        snake_command += ["--cores", threads]

    # add software deployment args
    if sdm:
        if sdm == "conda":
            snake_command += ["--sdm conda"]
            if prefix:
                snake_command += ["--conda-prefix", prefix]
        if sdm == "apptainer":
            snake_command += [
                f"--sdm apptainer --apptainer-args '-B {workflow_basedir}:{workflow_basedir}'"
            ]
            if prefix:
                snake_command += ["--apptainer-prefix", prefix]
    # add snakemake default args
    if snake_default:
        snake_command += snake_default

    # add any additional snakemake commands
    if snake_args:
        snake_command += list(snake_args)

    # allow double-handling of --profile
    if profile:
        snake_command += ["--profile", profile]

    # allow double-handling of --workflow-profile
    if workflow_profile:
        # copy system default if not present
        copy_config(
            os.path.join(workflow_profile, "config.yaml"),
            system_config=system_workflow_profile,
            log=log,
        )

        snake_command += ["--workflow-profile", workflow_profile]

    # Run Snakemake!!!
    snake_command = " ".join(str(s) for s in snake_command)
    msg_box("Snakemake command", errmsg=snake_command, log=log)
    if not subprocess.run(snake_command, shell=True).returncode == 0:
        msg("ERROR: Snakemake failed", log=log)
        sys.exit(1)
    else:
        msg("Snakemake finished successfully", log=log)
    return 0


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
            "--taxonkit",
            default=lambda: os.path.join(os.getcwd(), ".taxonkit"),
            help="Define path to taxonkit data-dir",
            type=click.Path(),
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
            "--sdm",
            type=click.Choice(["apptainer", "conda"]),
            default="apptainer",
            help="Software deplolyment method",
            show_default=True,
        ),
        click.option(
            "--prefix",
            default=None,
            help="Custom softawre deployment directory",
            callback=get_deployment_prefix,
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
    """Read simulation command line args
    include simulate args with the sim_options decorator.
    """
    options = [
        click.option(
            "--skip-shuffle",
            help="Skip fastq shuffling",
            type=bool,
            default=False,
            required=False,
        ),
        click.option(
            "--rotate",
            help="Number of times to shuffle genome start for circular assemblies (2 or more for circular)",
            type=int,
            default=1,
            show_default=True,
        ),
        click.option(
            "--tech",
            help="Sequencing technology",
            type=click.Choice(["illumina", "pacbio", "nanopore"], case_sensitive=False),
            default="illumina",
            show_default=True,
        ),
        click.option(
            "--paired",
            help="Illumina reads pairing",
            type=bool,
            default=True,
            show_default=True,
        ),
        click.option(
            "--tool",
            help="Which simulator to use",
            type=click.Choice(
                ["art", "pbsim3"],
                case_sensitive=False,
            ),
            default="art",
            show_default=True,
        ),
        click.option(
            "--error",
            help="Simulator error profile",
            type=str,
            default="HS25",
            show_default=True,
        ),
        click.option(
            "--custom-err",
            help="Path to art custom error profile basename",
            type=str,
            default=None,
        ),
        click.option(
            "--replicates",
            help="Number of replicates per sample",
            type=int,
            default=1,
            show_default=True,
        ),
        click.option(
            "--rep-sd",
            help="Standard deviation between replicate coverages (0: no deviation)",
            type=int,
            default=0,
            show_default=True,
        ),
        click.option(
            "--dist",
            help="Sample taxonomic abundances from a distribution",
            type=click.Choice(["none", "lognormal", "even"]),
            default="none",
            show_default=True,
        ),
        click.option(
            "--mu",
            help="Mu of lognormal dist",
            type=float,
            default=1,
        ),
        click.option(
            "--sigma",
            help="Sigma of lognormal dist",
            type=float,
            default=0,
        ),
        click.option(
            "--bases",
            help="Per sample base counts (ignored when coverage is set)",
            type=str,
            default="1G",
            show_default=True,
        ),
        click.option(
            "--model-path",
            help="Path to pbsim3 sequencing error models",
            type=click.Path(),
            default=snake_base(os.path.join("data", "models", "pbsim3")),
            show_default=True,
        ),
        click.option(
            "--model",
            help="Path to pbsim3 sequencing error models",
            type=click.Choice(["QSHMM-ONT-HQ", "QSHMM-ONT", "QSHMM-RSII"]),
        ),
        click.option(
            "--ratio",
            help="PBSIM3 substitution, insertion and deletion ratio",
            type=str,
            default="6:55:39",
            show_default=True,
        ),
        click.option(
            "--accuracy",
            help="Mean accuracy for long read sequencing",
            type=float,
            default=0.99,
            show_default=True,
        ),
        click.option(
            "--passes",
            help="Nmber of passes for pacbio multipass sequencing (2 or more for hifi)",
            type=int,
            default=1,
            show_default=True,
        ),
        click.option(
            "--min-len",
            help="Minimum read length for long read sequencing",
            type=int,
            default=100,
            show_default=True,
        ),
        click.option(
            "--max-len",
            help="Maximum read length for long read sequencing",
            type=int,
            default=1000000,
            show_default=True,
        ),
        click.option(
            "--sd-len",
            help="Standard read length deviation for long read sequencing",
            type=int,
            default=7000,
            show_default=True,
        ),
        click.option(
            "--mean-len",
            help="Mean read length for long and short read sequencing",
            type=int,
            default=150,
        ),
        click.option(
            "--frag-len",
            help="Fragment length for paired short read sequencing",
            type=int,
            default=200,
            show_default=True,
        ),
        click.option(
            "--frag-sd",
            help="Fragment length standard deviation for paired short read sequencing",
            type=int,
            default=10,
            show_default=True,
        ),
        click.option(
            "--bam/--no-bam",
            help="Generate gold standard bam files",
            default=False,
            show_default=True,
        ),
        click.option(
            "--ranks",
            help="Ranks to show in the taxonomic profile",
            type=str,
            default="superkingdom,phylum,class,order,family,genus,species,strain",
            show_default=True,
        ),
        click.option(
            "--seed",
            help="Seed for pseudo random number generator",
            type=int,
            default=1,
            show_default=True,
        ),
        click.option(
            "--errfree",
            help="Generate a zero sequencing errors SAM file",
            is_flag=True,
            default=True,
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
            "--limit",
            help="Limit number of genomes per query",
            type=str,
            default=1,
        ),
        click.option("--api-key", type=str, help="NCBI api-key", default=None),
        click.option(
            "--compressed",
            type=bool,
            help="Download compressed files",
            default=True,
            show_default=True,
        ),
        click.option(
            "--include",
            type=str,
            help="Comma seperated files to download : genome,rna,protein,cds,gff3,gtf,gbff,seq-report,none",
            default="genome,seq-report",
            show_default=True,
        ),
        click.option(
            "--source",
            type=click.Choice(["refseq", "genbank", "all"], case_sensitive=False),
            help="Download from refseq or genbank or both",
            default="all",
            show_default=True,
        ),
        click.option(
            "--taxon/--accession",
            help="Are queries taxa names or accession",
            type=bool,
            default=True,
            show_default=True,
        ),
        click.option(
            "--reference",
            type=bool,
            help="Limit to reference and representative genomes",
            default=True,
            show_default=True,
        ),
        click.option(
            "--assembly-level",
            help="Comma seperated list of assembly level: complete,chromosome,scaffold,contig",
            default=None,
            show_default=True,
        ),
        click.option(
            "--annotated",
            type=bool,
            help="Select annotated genomes only",
            default=False,
            show_default=True,
        ),
        click.option(
            "--atypical",
            type=bool,
            help="Exclude atypical genomes",
            default=True,
            show_default=True,
        ),
        click.option(
            "--mag",
            type=click.Choice(["exclude", "all", "only"], case_sensitive=False),
            help="Exclude, include or limit to metagenome assembled genomes",
            default="all",
            show_default=True,
        ),
        click.option(
            "--rank",
            help="taxonomic rank to filter by assemblies ",
            default=None,
            type=click.Choice(
                [
                    "superkingdom",
                    "phylum",
                    "class",
                    "order",
                    "family",
                    "genus",
                    "species",
                ],
                case_sensitive=False,
            ),
            show_default=True,
        ),
        click.option(
            "--nrank",
            help="Number of genomes per taxonomic rank",
            type=int,
            default=None,
            show_default=True,
        ),
    ]

    for option in reversed(options):
        func = option(func)
    return func
