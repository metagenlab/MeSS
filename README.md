# MeSS

[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)
[![Documentation Status](https://readthedocs.org/projects/metagenomic-sequence-simulator-mess/badge/?version=latest)](https://metagenomic-sequence-simulator-mess.readthedocs.io/en/latest/?badge=latest)
[![version](https://img.shields.io/conda/v/bioconda/mess?label=version&color=blue)](http://bioconda.github.io/recipes/mess/README.html)
[![downloads](https://img.shields.io/conda/dn/bioconda/mess.svg)](https://anaconda.org/bioconda/mess)

Metagenomic Sequence Simulator (MeSS) is a Snakemake pipeline, implemented using [Snaketool](https://github.com/beardymcjohnface/Snaketool), for simulating illumina, Oxford Nanopore (ONT) and Pacific Bioscience (PacBio) shotgun metagenomic samples.

For more details check [metagenomic-sequence-simulator-mess.readthedocs.io](https://metagenomic-sequence-simulator-mess.readthedocs.io)

## Installation

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/mess/README.html)

```sh
mamba create -n mess mess
```

## Usage

```sh
mess run -i <input> -o <outdir>
```

### Run minimal example

```sh
mess test
```

## Commands

```sh
        ___  ___     _____ _____
        |  \/  |    /  ___/  ___|
        | .  . | ___\ `--.\ `--.
        | |\/| |/ _ \`--. \`--. \
        | |  | |  __/\__/ /\__/ /
        \_|  |_/\___\____/\____/

Usage: mess [OPTIONS] COMMAND [ARGS]...

  For more options, run: mess command --help

Options:
  -v, --version  Show the version and exit.
  -h, --help     Show this message and exit.

Commands:
  run           Run MeSS
  download      Download assemblies
  simulate      Simulate reads from local fastas
  test          Run mess on test data
  hmp-template  Download and simulate healthy human microbiome templates
  config        Copy the system default config file
  citation      Print the citation(s) for this tool
```
