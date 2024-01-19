# Installation

The recommended way to install MeSS is with [mamba](https://github.com/mamba-org/mamba)

## Mamba

### Install mamba

```sh
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

### Configure channels for bioconda

```sh
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

### Install MeSS with mamba

```sh
mamba create -n mess mess
```

## Cloning the repo

`MeSS` cli can be installed via `pip`.

```sh
git clone https://github.com/metagenlab/MeSS
cd MeSS
pip install -e .
```
