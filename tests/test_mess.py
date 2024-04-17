import subprocess
from pathlib import Path
import pytest
import shutil
import os


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


test_data_path = Path("mess/test_data")
outdir = Path("test_out")
threads = 2


def remove_directory(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0
    Parameters
    ----------
    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None


def test_mess_cli():
    exec_command("mess -h")
    exec_command("mess run -h")
    exec_command("mess download -h")
    exec_command("mess simulate -h")
    exec_command("mess hmp-template -h")
    exec_command("mess config -h")
    exec_command("mess test -h")
    exec_command("mess citation")


def test_run():
    """run mess test"""
    cmd = f"mess test --threads {threads} --output {outdir}"
    exec_command(cmd)
    remove_directory(outdir)


cmd = (
    f"mess simulate",
    f"--threads {threads}",
    f"--input {test_data_path}/simulate_test.tsv",
    f"--fasta {test_data_path}/fastas",
    f"--output {outdir}",
)


def test_simulate_illumina():
    """mess simulate illumina reads"""
    exec_command(" ".join(cmd) + " --tech illumina --bam")
    remove_directory(outdir)


def test_simulate_nanopore():
    """mess simulate nanopore reads"""
    exec_command(" ".join(cmd) + " --tech nanopore --bam")
    remove_directory(outdir)


def test_simulate_pacbio():
    """mess simulate pacbio reads"""
    exec_command(" ".join(cmd) + " --tech pacbio --bam")
    remove_directory(outdir)