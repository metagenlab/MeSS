"""
Entrypoint for MeSS

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import rich_click as click

from .util import (
    copy_config,
    snake_base,
    run_snakemake,
    get_version,
    print_citation,
    common_options,
    download_options,
    sim_options,
)

click.rich_click.STYLE_COMMANDS_TABLE_COLUMN_WIDTH_RATIO = (1, 2)
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True

cmd_common_options = (
    {
        "name": "Common options",
        "options": ["--input", "--output", "--threads", "--taxonkit", "--help"],
    },
)
cmd_dl_options = (
    {
        "name": "Download options",
        "options": [
            "--api-key",
            "--limit",
            "--compressed",
            "--include",
            "--source",
            "--taxon",
            "--reference",
            "--assembly-level",
            "--annotated",
            "--atypical",
            "--mag",
            "--rank",
            "--nrank",
        ],
    },
)
cmd_sim_options = (
    {
        "name": "Common simulators options",
        "options": [
            "--tech",
            "--bases",
            "--tool",
            "--error",
            "--mean-len",
            "--bam",
            "--seed",
        ],
    },
    {
        "name": "art_illumina options",
        "options": ["--custom-err", "--paired", "--frag-len", "--frag-sd", "--errfree"],
    },
    {
        "name": "pbsim3 options",
        "options": [
            "--model",
            "--model-path",
            "--ratio",
            "--accuracy",
            "--passes",
            "--min-len",
            "--max-len",
            "--sd-len",
        ],
    },
    {
        "name": "Replicates options",
        "options": [
            "--replicates",
            "--rep-sd",
        ],
    },
    {
        "name": "Abundance distribution options",
        "options": [
            "--dist",
            "--mu",
            "--sigma",
        ],
    },
    {
        "name": "Taxonomic profile options",
        "options": [
            "--abundance",
            "--ranks",
        ],
    },
)
local_sim_options = (
    {
        "name": "Local genomes options",
        "options": [
            "--asm-summary",
            "--fasta",
            "--rotate",
        ],
    },
)

cmd_hmp_options = ({"name": "Hmp-template options", "options": ["--site", "--sample"]},)


skip_options = (
    {
        "name": "Skip options",
        "options": [
            "--skip-shuffle",
        ],
    },
)

snakemake_options = (
    {
        "name": "Snakemake options",
        "options": [
            "--configfile",
            "--profile",
            "--sdm",
            "--prefix",
            "--snake-default",
        ],
    },
)

mess_run = [
    options
    for options in cmd_common_options
    + cmd_dl_options
    + cmd_sim_options
    + skip_options
    + snakemake_options
]

mess_download = [
    options for options in cmd_common_options + cmd_dl_options + snakemake_options
]

mess_simulate = [
    options
    for options in cmd_common_options
    + local_sim_options
    + cmd_sim_options
    + skip_options
    + snakemake_options
]

hmp_options = [
    options
    for options in cmd_hmp_options
    + cmd_common_options
    + cmd_dl_options
    + cmd_sim_options
    + skip_options
    + snakemake_options
]


click.rich_click.OPTION_GROUPS = {
    "mess": [
        {
            "name": "Help",
            "options": ["--help", "--version"],
        }
    ],
    "mess run": mess_run,
    "mess download": mess_download,
    "mess simulate": mess_simulate,
    "mess hmp-template": hmp_options,
}

click.rich_click.COMMAND_GROUPS = {
    "mess": [
        {
            "name": "Main usage",
            "commands": ["run", "download", "simulate", "hmp-template"],
        },
        {
            "name": "Config/Test",
            "commands": [
                "config",
                "test",
            ],
        },
        {
            "name": "Cite",
            "commands": [
                "citation",
            ],
        },
    ]
}

version = get_version()


@click.group(
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.version_option(version, "-v", "--version", is_flag=True)
def cli():
    """
    \b
    ___  ___     _____ _____ 
    |  \\/  |    /  ___/  ___|
    | .  . | ___\\ `--.\\ `--. 
    | |\\/| |/ _ \\`--. \\`--. \\
    | |  | |  __/\\__/ /\\__/ /
    \\_|  |_/\\___\\____/\\____/
    \b
    For more options: 

    mess [command] --help
    """
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
Use conda:          mess run ... --sdm conda
Change defaults:    mess run ... --snake-default="-k --nolock"
Add Snakemake args: mess run ... --dry-run --keep-going --touch
"""

help_download = """
\n
Downloads genomes from taxons/accessions in input table(s)
\n
EXAMPLES:
mess download -i [input] -o [output]
"""

help_simulate = """
\n
Simulate reads from local fasta
\n
EXAMPLES:
mess simulate -i [input] --tech [tech] --fasta [fasta-dir] -o [output]
\n
or
\n 
mess simulate -i [input] --tech [tech] --asm-summary [asm-summary] -o [output]
"""
help_test = """
\n
Test command to run mess from assembly download to read simulation
\n
EXAMPLES:
mess test -o [output]
"""
help_hmp_templates = """
\n
Command to download and simulate healthy human microbiome templates
\n
EXAMPLES:
\n
mess hmp-template --site gut -o gut
\n
or
\n
mess hmp-template --site buccal_mucosa --sample SRS013506 -o SRS013506
"""


@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option(
    "-i",
    "--input",
    help="Path to input table(s)",
    type=click.Path(
        exists=True, file_okay=True, dir_okay=True, readable=True, resolve_path=True
    ),
    required=True,
)
@download_options
@common_options
@sim_options
def run(
    **kwargs,
):
    """Run MeSS workflow: download and simulate commands"""
    # Config to add or update in configfile
    merge_config = {"args": kwargs}

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
@click.option(
    "-i",
    "--input",
    help="Path to input table(s)",
    type=click.Path(
        exists=True, file_okay=True, dir_okay=True, readable=True, resolve_path=True
    ),
    required=True,
)
@download_options
@common_options
def download(
    **kwargs,
):
    """Download genomes with assembly_finder"""
    # Config to add or update in configfile
    merge_config = {"args": kwargs}

    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "download.smk")),
        merge_config=merge_config,
        **kwargs,
    )


@click.command(
    epilog=help_simulate,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option(
    "-i",
    "--input",
    help="Path to input table(s)",
    type=click.Path(
        exists=True, file_okay=True, dir_okay=True, readable=True, resolve_path=True
    ),
    required=True,
)
@sim_options
@click.option(
    "--asm-summary",
    help="Summary table with genome sizes, fasta paths, contig counts...",
    type=str,
)
@click.option(
    "--fasta",
    help="Path to local fasta directory if no path is set in summary or input table",
    type=str,
)
@common_options
def simulate(
    **kwargs,
):
    """Simulate reads from local fastas"""
    if kwargs["tech"] != "illumina":
        kwargs["mean_len"] = 9000
        if kwargs["tech"] == "pacbio":
            kwargs["ratio"] = "22:45:33"
            kwargs["model"] = "QSHMM-RSII"
            kwargs["accuracy"] = 0.85
            if kwargs["error"] == "hifi":
                kwargs["accuracy"] = 0.999
                kwargs["passes"] = 10
        if kwargs["tech"] == "nanopore":
            kwargs["ratio"] = "39:24:36"
            kwargs["accuracy"] = 0.99
            kwargs["model"] = "QSHMM-ONT-HQ"

    # Config to add or update in configfile
    merge_config = {"args": kwargs}
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "simulate.smk")),
        merge_config=merge_config,
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
    if kwargs["tech"] != "illumina":
        kwargs["mean_len"] = 9000
        if kwargs["tech"] == "pacbio":
            kwargs["ratio"] = "22:45:33"
            kwargs["model"] = "QSHMM-RSII"
            kwargs["accuracy"] = 0.85
            if kwargs["error"] == "hifi":
                kwargs["accuracy"] = 0.999
                kwargs["passes"] = 10
        if kwargs["tech"] == "nanopore":
            kwargs["ratio"] = "39:24:36"
            kwargs["accuracy"] = 0.99
            kwargs["model"] = "QSHMM-ONT-HQ"

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
    help="Choose microbiome site",
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
    if sample is None:
        kwargs["input"] = snake_base(os.path.join("data", "hmp_templates", site))
    else:
        kwargs["input"] = snake_base(
            os.path.join("data", "hmp_templates", site, f"{sample}.tsv")
        )
    kwargs["taxon"] = False
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
    cli()


if __name__ == "__main__":
    main()
