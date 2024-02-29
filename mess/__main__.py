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
)
from .util import (
    snake_base,
    get_version,
    print_citation,
    common_options,
    download_options,
    sim_options,
)


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

help_download = """
\b
downloads genomes from entries in input table(s)
mess download ...
\b
EXAMPLES:
mess download -i [input] -o [output]
"""

help_simulate = """
\b
simulate reads from local fasta
mess simulate ...
\b
EXAMPLES:
mess simulate -i [input] -o [output] --asm-summary [asm_summary] --tech [tech] 
"""
help_test = """
\b
test command to run mess from assembly download to read simulation
\b
EXAMPLES:
mess test -o [output]
"""
help_hmp_templates = """
\b
command to download and simulate healthy human microbiome templates
\b
EXAMPLES:
mess hmp_template --site gut -o gut
mess hmp_template --site buccal_mucosa --sample SRS013506 -o SRS013506
"""


@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("-i", "--input", help="Input sample sheet(s)", type=str, required=True)
@download_options
@common_options
@sim_options
def run(
    input,
    output,
    log,
    ncbi_email,
    ncbi_key,
    db,
    uid,
    rc,
    al,
    an,
    rank,
    nr,
    ete_db,
    nb,
    exclude,
    tech,
    tool,
    error,
    custom_err,
    replicates,
    rep_sd,
    dist,
    mu,
    sigma,
    bases,
    model,
    ratio,
    accuracy,
    passes,
    min_len,
    max_len,
    mean_len,
    sd_len,
    frag_len,
    frag_sd,
    bam,
    paired,
    seed,
    **kwargs,
):
    """Run MeSS"""
    # Config to add or update in configfile
    merge_config = {
        "args": {
            "input": input,
            "output": output,
            "log": log,
            "ncbi_email": ncbi_email,
            "ncbi_key": ncbi_key,
            "db": db,
            "uid": uid,
            "rc": rc,
            "al": al,
            "an": an,
            "rank": rank,
            "nr": nr,
            "ete_db": ete_db,
            "nb": nb,
            "exclude": exclude,
            "tech": tech,
            "tool": tool,
            "error": error,
            "custom_err": custom_err,
            "replicates": replicates,
            "rep_sd": rep_sd,
            "dist": dist,
            "mu": mu,
            "sigma": sigma,
            "bases": bases,
            "model": model,
            "ratio": ratio,
            "accuracy": accuracy,
            "passes": passes,
            "min_len": min_len,
            "max_len": max_len,
            "mean_len": mean_len,
            "sd_len": sd_len,
            "frag_len": frag_len,
            "frag_sd": frag_sd,
            "bam": bam,
            "paired": paired,
            "seed": seed,
        }
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "Snakefile")),
        merge_config=merge_config,
        **kwargs,
    )


@click.command(
    epilog=help_download,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("-i", "--input", help="Input file/directory", type=str, required=True)
@download_options
@common_options
def download(
    input,
    output,
    log,
    ncbi_email,
    ncbi_key,
    db,
    uid,
    rc,
    al,
    an,
    rank,
    nr,
    ete_db,
    nb,
    exclude,
    **kwargs,
):
    """Download assemblies"""
    # Config to add or update in configfile
    merge_config = {
        "args": {
            "input": input,
            "output": output,
            "log": log,
            "ncbi_email": ncbi_email,
            "ncbi_key": ncbi_key,
            "db": db,
            "uid": uid,
            "rc": rc,
            "al": al,
            "an": an,
            "rank": rank,
            "nr": nr,
            "ete_db": ete_db,
            "nb": nb,
            "exclude": exclude,
        }
    }
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "download.smk")),
        merge_config=merge_config,
        log=log,
        **kwargs,
    )


@click.command(
    epilog=help_simulate,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("-i", "--input", help="path to sample sheet(s)", type=str, required=True)
@sim_options
@click.option(
    "--skip-fa-proc",
    help="skip fasta processing (decompressing, renaming headers, merge/split contigs)",
    type=bool,
    default=False,
    required=False,
)
@click.option(
    "--asm-summary",
    help="path to assembly summary table (contains taxid, genome_size, path to fasta)",
    type=str,
    required=True,
)
@common_options
def simulate(
    input,
    output,
    log,
    skip_fa_proc,
    skip_shuffle,
    compressed,
    asm_summary,
    tech,
    tool,
    error,
    model,
    ratio,
    custom_err,
    replicates,
    rep_sd,
    dist,
    mu,
    sigma,
    bases,
    accuracy,
    passes,
    min_len,
    max_len,
    sd_len,
    mean_len,
    frag_len,
    frag_sd,
    bam,
    paired,
    seed,
    **kwargs,
):
    if tech != "illumina" and error:
        if tech == "pacbio":
            ratio = "22:45:33"
            model = "QSHMM-RSII"
            if error == "hifi":
                accuracy = 0.99
                passes = 10

        if tech == "nanopore":
            ratio = "39:24:36"
            if error == "r10.4":
                accuracy = 0.99
                model = "QSHMM-ONT-HQ"
            if error == "r10.3":
                accuracy = 0.95
                model = "QSHMM-ONT"

    """Simulate reads from local fastas"""
    # Config to add or update in configfile
    merge_config = {
        "args": {
            "input": input,
            "output": output,
            "log": log,
            "skip_fa_proc": skip_fa_proc,
            "skip_shuffle": skip_shuffle,
            "compressed": compressed,
            "tech": tech,
            "tool": tool,
            "error": error,
            "model": model,
            "ratio": ratio,
            "custom_err": custom_err,
            "asm_summary": asm_summary,
            "replicates": replicates,
            "rep_sd": rep_sd,
            "dist": dist,
            "mu": mu,
            "sigma": sigma,
            "bases": bases,
            "accuracy": accuracy,
            "passes": passes,
            "min_len": min_len,
            "max_len": max_len,
            "sd_len": sd_len,
            "mean_len": mean_len,
            "frag_len": frag_len,
            "frag_sd": frag_sd,
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


@click.command(
    epilog=help_test,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@download_options
@sim_options
@common_options
def test(**kwargs):
    """Run mess on test data"""
    # Config to add or update in configfile
    kwargs["input"] = snake_base(os.path.join("test_data", "minimal_test.tsv"))
    merge_config = {"args": kwargs}
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "Snakefile")),
        merge_config=merge_config,
        **kwargs,
    )


@click.command(
    epilog=help_hmp_templates,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option(
    "--site",
    help="choose microbiome site",
    type=click.Choice(["gut", "buccal_mucosa", "vagina", "throat"]),
    required=True,
    default="gut",
)
@click.option(
    "--sample", help="choose sample template", type=str, required=False, default=None
)
@download_options
@sim_options
@common_options
def hmp_template(site, sample, **kwargs):
    """Download and simulate healthy human microbiome templates"""
    # Config to add or update in configfile
    if sample == None:
        kwargs["input"] = snake_base(os.path.join("data", site))
    else:
        kwargs["input"] = snake_base(os.path.join("data", site, f"{sample}.tsv"))
    merge_config = {"args": kwargs}
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "Snakefile")),
        merge_config=merge_config,
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
cli.add_command(test)
cli.add_command(hmp_template)
cli.add_command(config)
cli.add_command(citation)


def main():
    print_splash()
    cli()


if __name__ == "__main__":
    main()
