#Parts of this wrapper are inspired by atlas, an easy-to-use metagenomic pipeline based on snakemake.
#Go check it out on https://github.com/metagenome-atlas/atlas
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
    short_help="run all MeSS pipeline steps"
)
@click.option("-f",
    "--config-file",
    type=click.Path(exists=True,resolve_path=True),
    help="path to MeSS configfile",
)
@click.option("-p",
    "--conda-prefix",
    type=click.Path(exists=True,resolve_path=True),
    help="path to conda environment",
)
@click.option(
    "-n",
    "--dryrun_status",
    is_flag=True,
    default=False,
    show_default=True,
    help="Snakemake dryrun to see the scheduling plan",
)
@click.option("-nr",
    "--ncbi_requests",
    type=int,
    help="number of ncbi requests, check https://github.com/metagenlab/mess for more details",
    default=3,
)
@click.option("-ns",
    "--nb_sim",
    type=int,
    help="number of genomes to simulate in parallel, check https://github.com/metagenlab/mess for more details",
    default=2,
)
@click.option("-nc",
    "--nb_cat",
    type=int,
    help="number of genomes to concatenate in parallel, check https://github.com/metagenlab/mess for more details",
    default=2,
)
@click.option("-c",
    "--cores",
    type=int,
    help="number of cores to allow for the workflow",
    default=5,
)
@click.option("-st",
    "--seq-tech",
    help="sequencing technologie for simulation",
)

@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def run_workflow(conda_prefix,config_file,dryrun_status,ncbi_requests,nb_sim,nb_cat,cores,snakemake_args,seq_tech):
    """
    Runs MeSS pipeline with all steps
    """
    if config_file is None:#if config-file path is not set, assume that config-file is present where the command is ran
        config_file = 'config.yml'

    if not os.path.exists(config_file):
        logging.critical(f"config-file not found: {config_file}\n"
                         "make sure your config.yml is in your working directory")
        sys.exit(1)
    if dryrun_status:
        dryrun='-n'
    else:
        dryrun=''
    if conda_prefix is None:
        conda_prefix=os.environ['CONDA_PREFIX']
    if not os.path.exists(conda_prefix):
        logging.critical(f"conda env path not found: {config_file}\n")
        sys.exit(1)
    if seq_tech=='illumina':
        cmd = (
            f"snakemake --snakefile {get_snakefile()} --configfile {config_file} --use-conda --conda-prefix {conda_prefix} "
            f" --resources ncbi_requests={ncbi_requests} nb_simulation={nb_sim} parallel_cat={nb_cat} --cores {cores} all_sim {dryrun} {' '.join(snakemake_args)}")
    if seq_tech == 'longreads':
        cmd = (
            f"snakemake --snakefile {get_snakefile()} --configfile {config_file} --use-conda --conda-prefix {conda_prefix} --use-singularity "
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
