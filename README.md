# Metagenomic Sequence Simulator (MeSS)

[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)
[![license](https://img.shields.io/github/license/metagenlab/mess.svg)](https://github.com/metagenlab/MeSS/blob/main/LICENSE)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/mess/README.html)
[![version](https://img.shields.io/conda/vn/bioconda/mess?color=blue)](http://bioconda.github.io/recipes/mess/README.html)
[![downloads](https://img.shields.io/conda/dn/bioconda/mess.svg)](https://anaconda.org/bioconda/mess)

[![tests](https://github.com/metagenlab/MeSS/actions/workflows/unit-tests.yml/badge.svg)](https://github.com/metagenlab/MeSS/actions/workflows/unit-tests.yml)
[![docs](https://github.com/metagenlab/MeSS/actions/workflows/build-docs.yml/badge.svg)](https://github.com/metagenlab/MeSS/actions/workflows/build-docs.yml)
[![docker](https://github.com/metagenlab/MeSS/actions/workflows/docker-publish.yml/badge.svg)](https://github.com/metagenlab/MeSS/actions/workflows/docker-publish.yml)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13365501.svg)](https://zenodo.org/doi/10.5281/zenodo.13365501)

The Metagenomic Sequence Simulator (MeSS) is a [Snakemake](https://github.com/snakemake/snakemake) pipeline, implemented using [Snaketool](https://github.com/beardymcjohnface/Snaketool), for simulating illumina, Oxford Nanopore (ONT) and Pacific Bioscience (PacBio) shotgun metagenomic samples.

## :mag: Overview

MeSS takes as input NCBI taxa or local genome assemblies to generate either long (PacBio or ONT) or short (illumina) reads. In addition to reads, MeSS optionally generates bam alignment files and taxonomic + sequence abundances in [CAMI format](https://github.com/bioboxes/rfc/blob/master/data-format/profiling.mkd).

```mermaid
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
style genome_download color:#15161a

input --> distchoice
subgraph community_design["`**community design**`"]
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
style community_design color:#15161a
style community_design color:#15161a

fasta --> simulator
depth --> simulator

simulator["read simulator
(art_illumina, pbsim3...)"]
simulator --> bam
simulator --> fastq
simulator --> CAMI-profile

%% subgraph color fills
classDef red fill:#faeaea,color:#fff,stroke:#333;
classDef blue fill:#eaecfa,color:#fff,stroke:#333;
class genome_download blue

class community_design red
```

## :books: Documentation

More details can be found in the [documentation](https://metagenlab.github.io/MeSS/)

## :zap: Quick start

### :gear: Installation

- Conda ([Miniforge](https://github.com/conda-forge/miniforge))

```sh
conda create -n mess mess
```

- Docker

```sh
docker pull ghcr.io/metagenlab/mess:latest
```

- From source

```sh
git clone https://github.com/metagenlab/MeSS.git
pip install -e MeSS
```

### :page_facing_up: Usage

#### :arrow_right: Input

Let's simulate two metagenomic samples with the following taxa and read counts in `samples.tsv`:
| sample | taxon | reads |
| --- | --- | --- |
| sample1 | 487 | 174840 |
| sample1 | 727 | 90679 |
| sample1 | 729 | 13129 |
| sample2 | 28132 | 147863 |
| sample2 | 199 | 147545 |
| sample2 | 729 | 131300 |

#### :rocket: Command

```sh
mess run -i samples.tsv
```

> [!IMPORTANT] 
> [Apptainer](https://apptainer.org/) is the default and recommended dependency deployment method for maximum reproducibility ! If you would like to use conda you can specify `--sdm conda`.

#### :card_index_dividers: Outputs

- Downloaded genomes in `mess_out/assembly_finder/download`

```sh
┣ 📂GCF_000144405.1
┃ ┗ 📜GCF_000144405.1_ASM14440v1_genomic.fna.gz
┣ 📂GCF_001298465.1
┃ ┗ 📜GCF_001298465.1_ASM129846v1_genomic.fna.gz
┣ 📂GCF_016127215.1
┃ ┗ 📜GCF_016127215.1_ASM1612721v1_genomic.fna.gz
┣ 📂GCF_020736045.1
┃ ┗ 📜GCF_020736045.1_ASM2073604v1_genomic.fna.gz
┣ 📂GCF_022869645.1
 ┗ 📜GCF_022869645.1_ASM2286964v1_genomic.fna.gz
```

- Simulated reads in `mess_out/fastq`

```sh
┣ 📜sample1_R1.fq.gz
┣ 📜sample1_R2.fq.gz
┣ 📜sample2_R1.fq.gz
┗ 📜sample2_R2.fq.gz
```

> [!TIP]
> By default `mess` outputs paired illumina reads with the Hiseq25k error profile. Other outputs, and error profiles are described [here](https://metagenlab.github.io/MeSS/guide/output/) and [here](https://metagenlab.github.io/MeSS/tutorials/seqtech/)

#### :bar_chart: Resources usage

Using [`samples.tsv`](#arrow_right-input), `mess` runs in under 2min, while using around 1.8GB of physical RAM

| task_id | hash      | native_id | name     | status    | exit | submit                  | duration | realtime | %cpu   | peak_rss | peak_vmem | rchar  | wchar  |
| ------- | --------- | --------- | -------- | --------- | ---- | ----------------------- | -------- | -------- | ------ | -------- | --------- | ------ | ------ |
| 1       | fe/03c2bc | 62286     | MESS (1) | COMPLETED | 0    | 2024-09-04 12:41:15.820 | 1m 50s   | 1m 50s   | 111.5% | 1.8 GB   | 9 GB      | 3.5 GB | 2.4 GB |
| 1       | ff/0d03b1 | 73355     | MESS (1) | COMPLETED | 0    | 2024-09-04 12:55:12.903 | 1m 52s   | 1m 52s   | 112.6% | 1.7 GB   | 8.8 GB    | 3.5 GB | 2.4 GB |
| 1       | 07/d352bf | 83576     | MESS (1) | COMPLETED | 0    | 2024-09-04 12:57:30.600 | 1m 50s   | 1m 50s   | 113.2% | 1.7 GB   | 8.9 GB    | 3.5 GB | 2.4 GB |

> [!NOTE]
> Average resources usage measured 3 times with one CPU (using [nextflow](https://github.com/nextflow-io/nextflow), excluding dependency deployment time).

More details in the [resource usage documentation](https://metagenlab.github.io/MeSS/benchmarks/resource-usage/)

## :fire: Features

### :dna: Multi sequencing technology choice

- Illumina

```sh
mess test --tech illumina
```

- Nanopore

```sh
mess test --tech nanopore
```

- PacBio

```sh
mess test --tech pacbio
```

### :white_check_mark: BAMs and taxonomic profiles

```sh
mess test --bam
```

### :o: Circular genomes

```sh
mess test --rotate 3
```


## :sos: Help

All command-line options at described [here](https://metagenlab.github.io/MeSS/commands/)

![`mess -h`](docs/images/mess-help.svg)
