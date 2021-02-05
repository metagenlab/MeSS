import logging
import os
import sys
import click
import subprocess
from mess import __version__
logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M",
    format="[%(asctime)s %(levelname)s] %(message)s",
)

def get_snakefile():
    thisdir = os.path.abspath(os.path.dirname(__file__))
    sf = os.path.join(thisdir,'scripts','Snakefile')
    if not os.path.exists(sf):
        sys.exit(f"Unable to locate the Snakemake workflow file at {sf}")
    return sf

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """
    MeSS is a snakemake workflow used for simulating metagenomic mock communities
    """
@cli.command(
    "run",
    context_settings=dict(ignore_unknown_options=True),
    short_help="run MeSS simulator"
)
@click.option("-w",
    "--working-dir",
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    help="location to MeSS.",
    default="."
)
@click.option("-c",
    "--config-file",
    type=click.Path(exists=True,resolve_path=True),
    help="path to MeSS configfile",)

@click.option(
    "-n",
    "--dryrun_status",
    is_flag=True,
    default=False,
    show_default=True,
    help="Snakemake dryrun",
)
@click.option("-r",
    "--ncbi_requests",
    type=int,
    help="number of ncbi requests, check README.md for more details",
    default=3,
)
@click.option("-s",
    "--nb_sim",
    type=int,
    help="number of genomes to simulate in parallel, check README.md for more details",
    default=2,
)
@click.option("-k",
    "--nb_cat",
    type=int,
    help="number of genomes to concatenate in parallel, check README.md for more details",
    default=2,
)
@click.option("-c",
    "--cores",
    type=int,
    help="number of cores to allow for the workflow",
    default=5,
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def run_workflow(working_dir,config_file,dryrun_status,ncbi_requests,nb_sim,nb_cat,cores,snakemake_args):
    """
    Runs MeSS pipeline with all steps
    """
    if config_file is None:
        config_file = os.path.join(working_dir, 'config.yml')

    if not os.path.exists(config_file):
        logging.critical(f"config-file not found: {config_file}\n"
                         "make sure your config.yml is in your working directory")
        sys.exit(1)
    if dryrun_status:
        dryrun='-n'
    else:
        dryrun=''

    cmd = (
        f"snakemake --snakefile {get_snakefile()} --configfile {config_file} --use-conda --conda-prefix ${{CONDA_DEFAULT_ENV}} "
        f" --resources ncbi_requests={ncbi_requests} nb_simulation={nb_sim} parallel_cat={nb_cat} --cores {cores} all_sim {dryrun} {' '.join(snakemake_args)}")

    logging.info("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
        exit(1)
if __name__ == "__main__":
    cli()
