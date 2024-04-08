"""
Entrypoint for MeSS

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import rich_click as click


from snaketool_utils.cli_utils import (
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

click.rich_click.STYLE_COMMANDS_TABLE_COLUMN_WIDTH_RATIO = (1, 2)
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True

mess_common_options = (
    {
        "name": "Common options",
        "options": ["--input", "--output", "--threads", "--help"],
    },
)
mess_download_options = (
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
mess_simulate_options = (
    {
        "name": "Simulators options",
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
        "name": "ART illumina options",
        "options": ["--custom-err", "--paired", "--frag-len", "--frag-sd"],
    },
    {
        "name": "PBSIM3 options",
        "options": [
            "--model",
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
)
mess_local_sim_options = (
    {
        "name": "Mess simulate options",
        "options": ["--asm-summary"],
    },
)

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
            "--use-conda",
            "--conda-prefix",
            "--snake-default",
        ],
    },
)

mess_run = [
    options
    for options in mess_common_options
    + mess_download_options
    + mess_simulate_options
    + skip_options
    + snakemake_options
]

mess_download = [
    options
    for options in mess_common_options + mess_download_options + snakemake_options
]

mess_simulate = [
    options
    for options in mess_common_options
    + mess_local_sim_options
    + mess_simulate_options
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
    "mess hmp-template": [
        {"name": "Hmp-template options", "options": ["--site", "--sample"]}
    ],
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
    |  \/  |    /  ___/  ___|
    | .  . | ___\ `--.\ `--. 
    | |\/| |/ _ \`--. \`--. \\
    | |  | |  __/\__/ /\__/ /
    \_|  |_/\___\____/\____/
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
@click.option("-i", "--input", help="Input file/directory", type=str, required=True)
@download_options
@common_options
def download(
    **kwargs,
):
    """Download genomes"""
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
@click.option("-i", "--input", help="path to sample sheet(s)", type=str, required=True)
@sim_options
@click.option(
    "--asm-summary",
    help="Path to assembly summary table (contains tax_id, genome_size, path to fasta)",
    type=str,
    required=True,
)
@common_options
def simulate(
    **kwargs,
):
    """Simulate fastq from local genomes"""
    if kwargs["tech"] != "illumina" and kwargs["error"]:
        if kwargs["tech"] == "pacbio":
            kwargs["ratio"] = "22:45:33"
            kwargs["model"] = "QSHMM-RSII"
            if kwargs["error"] == "hifi":
                kwargs["accuracy"] = 0.99
                kwargs["passes"] = 10

        if kwargs["tech"] == "nanopore":
            kwargs["ratio"] = "39:24:36"
            if kwargs["error"] == "r10.4":
                kwargs["accuracy"] = 0.99
                kwargs["model"] = "QSHMM-ONT-HQ"
            if kwargs["error"] == "r10.3":
                kwargs["accuracy"] = 0.95
                kwargs["model"] = "QSHMM-ONT"

    """Simulate reads from local fastas"""
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
    if kwargs["tech"] != "illumina" and kwargs["error"]:
        if kwargs["tech"] == "pacbio":
            kwargs["ratio"] = "22:45:33"
            kwargs["model"] = "QSHMM-RSII"
            if kwargs["error"] == "hifi":
                kwargs["accuracy"] = 0.99
                kwargs["passes"] = 2

        if kwargs["tech"] == "nanopore":
            kwargs["ratio"] = "39:24:36"
            if kwargs["error"] == "r10.4":
                kwargs["accuracy"] = 0.99
                kwargs["model"] = "QSHMM-ONT-HQ"
            if kwargs["error"] == "r10.3":
                kwargs["accuracy"] = 0.95
                kwargs["model"] = "QSHMM-ONT"

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
    cli()


if __name__ == "__main__":
    main()
