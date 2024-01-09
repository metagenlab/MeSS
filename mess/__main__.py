"""
Entrypoint for MeSS

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import click

from snaketool_utils.cli_utils import (
    OrderedCommands,
    run_snakemake,
    copy_config,
    echo_click,
)


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
    Define common command line args here, and include them with the @common_options decorator below.
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
        click.argument("snake_args", nargs=-1, type=click.UNPROCESSED),
    ]
    for option in reversed(options):
        func = option(func)
    return func


version = get_version()


def print_splash():
    click.echo(
        """
        ___  ___     _____ _____ 
        |  \/  |    /  ___/  ___|
        | .  . | ___\ `--.\ `--. 
        | |\/| |/ _ \`--. \`--. \\
        | |  | |  __/\__/ /\__/ /
        \_|  |_/\___\____/\____/
        """
    )


@click.group(
    cls=OrderedCommands, context_settings=dict(help_option_names=["-h", "--help"])
)
@click.version_option(version, "-v", "--version", is_flag=True)
def cli():
    """
    For more options, run:
    mess command --help"""
    pass


help_msg_extra = """
\b
CLUSTER EXECUTION:
mess run ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           mess run --input [file]
Specify threads:    mess run ... --threads [threads]
Disable conda:      mess run ... --no-use-conda 
Change defaults:    mess run ... --snake-default="-k --nolock"
Add Snakemake args: mess run ... --dry-run --keep-going --touch
Specify targets:    mess run ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""


@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("-i", "--input", help="Input file/directory", type=str, required=True)
@common_options
def run(**kwargs):
    """Run MeSS"""
    # Config to add or update in configfile
    merge_config = {"args": kwargs}

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "Snakefile")),
        merge_config=merge_config,
        **kwargs,
    )


@click.command()
@click.option("-i", "--input", help="Input file/directory", type=str, required=True)
@click.option("--ncbi_key", help="NCBI key", type=str, required=False)
@click.option("--ncbi_email", help="NCBI email", type=str, required=False)
@common_options
def download(input, output, log, ncbi_email, ncbi_key, **kwargs):
    """download assemblies"""
    # Config to add or update in configfile
    merge_config = {
        "args": {
            "input": input,
            "output": output,
            "log": log,
            "ncbi_email": ncbi_email,
            "ncbi_key": ncbi_key,
        }
    }
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "download.smk")),
        merge_config=merge_config,
        **kwargs,
    )


profiles = {
    "illumina": {
        "art": [
            "MSv3",
            "GA1",
            "GA2",
            "HS10",
            "HS20",
            "HS25",
            "HSXn",
            "HSXt",
            "MinS",
            "MSv1",
            "NS50",
            "custom",
        ]
    },
    "pacbio": {"pbsim3": ["clr", "hifi"]},
    "nanopore": {"pbsim3": ["r10.4", "r10.3"]},
}


@click.command()
@click.option("-i", "--input", help="path to sample sheet(s)", type=str, required=True)
@click.option(
    "--fasta",
    help="path to fasta file/directory",
    type=str,
    required=True,
)
@click.option(
    "--asm-summary",
    help="path to assembly summary table (contains taxid, genome_size and contig_counts)",
    type=str,
    required=True,
)
@click.option(
    "--tech",
    help="sequencing technology",
    type=click.Choice(["illumina", "pacbio", "nanopore"], case_sensitive=False),
    default="illumina",
    show_default=True,
)
@click.option(
    "--paired",
    help="illumina reads pairing",
    type=bool,
    default=True,
    show_default=True,
)
@click.option(
    "--tool",
    help="choose simulator",
    type=click.Choice(
        ["art", "pbsim3"],
        case_sensitive=False,
    ),
    default="art",
    show_default=True,
)
@click.option(
    "--err-profile",
    help="choose simulator profile",
    type=str,
    default="MSv3",
    show_default=True,
)
@click.option(
    "--replicates",
    help="number of replicates per sample",
    type=int,
    default=1,
    show_default=True,
)
@click.option(
    "--rep-sd",
    help="standard deviation between replicate coverages (0: no deviation)",
    type=int,
    default=0,
    show_default=True,
)
@click.option(
    "--dist",
    help="sample genome coverage values from a distribution",
    type=click.Choice(["none", "lognormal", "even"]),
    default="none",
    show_default=True,
)
@click.option(
    "--mu",
    help="mu of lognormal dist",
    type=int,
    default=1,
)
@click.option(
    "--sigma",
    help="sigma of lognormal dist",
    type=int,
    default=0,
)
@click.option(
    "--total-bases",
    help="total amount of bases to simulate (ignored when coverage is set)",
    type=str,
    default="1G",
    show_default=True,
)
@click.option(
    "--accuracy",
    help="mean accuracy for long read sequencing",
    type=float,
    default=0.9,
    show_default=True,
)
@click.option(
    "--passes",
    help="number of passes for pacbio multipass sequencing (2 or more for hifi)",
    type=int,
    default=1,
    show_default=True,
)
@click.option(
    "--min-len",
    help="minimum read length for long read sequencing",
    type=int,
    default=100,
    show_default=True,
)
@click.option(
    "--max-len",
    help="maximum read length for long read sequencing",
    type=int,
    default=300,
    show_default=True,
)
@click.option(
    "--mean-len",
    help="mean read length for long and short read sequencing",
    type=int,
    default=150,
    show_default=True,
)
@click.option(
    "--sd-len",
    help="standard read length deviation for long read sequencing",
    type=int,
    default=50,
    show_default=True,
)
@click.option(
    "--bam",
    help="generate gold standard alignment bam files",
    type=bool,
    default=False,
    show_default=True,
)
@click.option(
    "--seed",
    help="seed for pseudo random number generator",
    type=int,
    default=42,
    show_default=True,
)
@common_options
def simulate(
    input,
    output,
    log,
    fasta,
    asm_summary,
    tech,
    err_profile,
    replicates,
    rep_sd,
    dist,
    mu,
    sigma,
    total_bases,
    accuracy,
    passes,
    min_len,
    max_len,
    mean_len,
    sd_len,
    bam,
    paired,
    seed,
    **kwargs
):
    """simulate reads"""
    # Config to add or update in configfile
    merge_config = {
        "args": {
            "input": input,
            "output": output,
            "log": log,
            "fasta": fasta,
            "tech": tech,
            "err_profile": err_profile,
            "asm_summary": asm_summary,
            "replicates": replicates,
            "rep_sd": rep_sd,
            "dist": dist,
            "mu": mu,
            "sigma": sigma,
            "total_bases": total_bases,
            "accuracy": accuracy,
            "passes": passes,
            "min_len": min_len,
            "max_len": max_len,
            "mean_len": mean_len,
            "sd_len": sd_len,
            "bam": bam,
            "paired": paired,
            "seed": seed,
        }
    }
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "simulate.smk")),
        merge_config=merge_config,
        log=log,
        **kwargs,
    )


@click.command()
@common_options
def config(configfile, system_config, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile, system_config=system_config)


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


cli.add_command(run)
cli.add_command(download)
cli.add_command(simulate)
cli.add_command(config)
cli.add_command(citation)


def main():
    print_splash()
    cli()


if __name__ == "__main__":
    main()
