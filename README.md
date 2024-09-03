# Metagenomic Sequence Simulator (MeSS)

[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)
[![license](https://img.shields.io/github/license/metagenlab/mess.svg)](https://github.com/metagenlab/MeSS/blob/main/LICENSE)
[![version](https://img.shields.io/conda/vn/bioconda/mess?color=blue)](http://bioconda.github.io/recipes/mess/README.html)
[![downloads](https://img.shields.io/conda/dn/bioconda/mess.svg)](https://anaconda.org/bioconda/mess)

[![tests](https://github.com/metagenlab/MeSS/actions/workflows/unit-tests.yml/badge.svg)](https://github.com/metagenlab/MeSS/actions/workflows/unit-tests.yml)
[![docs](https://github.com/metagenlab/MeSS/actions/workflows/build-docs.yml/badge.svg)](https://github.com/metagenlab/MeSS/actions/workflows/build-docs.yml)
[![docker](https://github.com/metagenlab/MeSS/actions/workflows/docker-publish.yml/badge.svg)](https://github.com/metagenlab/MeSS/actions/workflows/docker-publish.yml)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13365501.svg)](https://zenodo.org/doi/10.5281/zenodo.13365501)


The Metagenomic Sequence Simulator (MeSS) is a [Snakemake](https://github.com/snakemake/snakemake) pipeline, implemented using [Snaketool](https://github.com/beardymcjohnface/Snaketool), for simulating illumina, Oxford Nanopore (ONT) and Pacific Bioscience (PacBio) shotgun metagenomic samples.

## :memo: Overview

MeSS takes as input NCBI taxa or local genome assemblies to generate either long (PacBio or ONT) or short (illumina) reads. In addition to reads, MeSS optionally generates bam alignment files and taxonomic + sequence abundances in [CAMI format](https://github.com/bioboxes/rfc/blob/master/data-format/profiling.mkd).

``` mermaid
%%{init: {'theme':'forest'}}%%
flowchart LR
input["samples.tsv 
or 
samples/*.tsv"] --> taxons

subgraph genome_download["genome download"]
dlchoice{download ?}
taxons["taxons or
accesions"] --> dlchoice
dlchoice -->|yes| assembly_finder
dlchoice -->|no| fasta 
assembly_finder --> fasta
end

input --> distchoice
subgraph community_design["community design"]
distchoice{draw distribution ?}
distchoice -->|yes| dist["distribution 
(lognormal, even)"]
dist --> abundances
distchoice -->|no| reads
distchoice -->|no| bases
distchoice -->|no| abundances
depth["coverage depth"]
reads --> depth
bases --> depth
abundances["abundances 
(sequence, taxonomic)"] --> depth 
end
fasta --> simulator
depth --> simulator

simulator["read simulator 
(art_illumina, pbsim3...)"]
simulator --> bam
simulator --> fastq
simulator --> CAMI-profile

%% colors
style genome_download color:black
style community_design color:black
classDef red fill:#faeaea,color:#fff,stroke:#333;
classDef blue fill:#eaecfa,color:#fff,stroke:#333;
class genome_download blue
class community_design red
```
## :books: Documentation 

More details can be found in the [documentation](https://metagenlab.github.io/MeSS/)

## :zap: Quick start 
### Installation

#### Mamba

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/mess/README.html)

```sh
mamba create -n mess mess
```

#### Docker

```sh
docker pull ghcr.io/metagenlab/mess:latest
```

#### From source 

```sh
git clone https://github.com/metagenlab/MeSS.git
pip install -e MeSS
```

### Usage


#### Download and simulate

Using the following file [minimal_test.tsv](https://github.com/metagenlab/MeSS/blob/main/mess/test_data/minimal_test.tsv)

```sh
mess run -i minimal_test.tsv 
```

#### Simulate from local fasta

Download the [fasta directory](https://github.com/metagenlab/MeSS/tree/main/mess/test_data/fastas) and [table](https://github.com/metagenlab/MeSS/blob/main/mess/test_data/simulate_test.tsv)

```sh
mess simulate -i simulate_test.tsv --fasta fasta 
```

## :sos: Help

More details on command-line options in the [doc](https://metagenlab.github.io/MeSS/commands/)

![`mess -h`](docs/images/mess-help.svg)
