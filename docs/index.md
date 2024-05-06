# Introduction

[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![version](https://img.shields.io/conda/v/bioconda/mess?label=version&color=blue)](http://bioconda.github.io/recipes/mess/README.html)
[![downloads](https://img.shields.io/conda/dn/bioconda/mess.svg)](https://anaconda.org/bioconda/mess)

The Metagenomic Sequence Simulator (MeSS) is a [Snakemake](https://github.com/snakemake/snakemake) pipeline, implemented using [Snaketool](https://github.com/beardymcjohnface/Snaketool), for simulating illumina, Oxford Nanopore (ONT) and Pacific Bioscience (PacBio) shotgun metagenomic samples.

<div class="grid cards" markdown>

-   :material-clock-fast:{ .lg .middle } __Set up in few minutes__

    ---

    Install [`mess`](https://github.com/metagenlab/MeSS) with [`mamba`](https://github.com/mamba-org/mamba) and get started within minutes

    [:octicons-arrow-right-24: Getting started](getting_started/index.md)

-   :material-swap-vertical-variant:{ .lg .middle } __Streamlined pipeline__

    ---

    Easily download genomes and simulate reads in a single command
    
    [:octicons-arrow-right-24: Step-by-step guide](guide/index.md)

-   :dna:{ .lg .middle } __Multi-sequencing technologies__

    ---

    Simulate sequencing runs for both long and short reads

    [:octicons-arrow-right-24: Tutorials](tutorials/index.md#sequencing-technologies)

-   :material-speedometer:{ .lg .middle } __Fast and scalable__

    ---

    MeSS can generate complex communities in a few minutes

    [:octicons-arrow-right-24: Benchmarks](benchmarks/index.md)



</div>

## Overview

MeSS takes as input NCBI taxa or local genome assemblies to generate either long (PacBio or ONT) or short (illumina) reads. In addition to reads, MeSS optionally generates bam alignment files and taxonomic profiles in [bioboxes format](https://github.com/bioboxes/rfc).

![overview](images/workflow.svg)

## Installation

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/mess/README.html)

```sh
mamba create -n mess mess
```

## Usage

![`mess -h`](images/mess-help.svg)
