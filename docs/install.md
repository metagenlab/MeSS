## Mamba <small>(recommended)</small>

Install mamba via miniforge

```sh
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

Next, configure strict channel priority:

```sh
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```
Install with mamba

```sh
mamba create -n mess mess
```

## Git
Install Mamba as shown before, clone the MeSS repo and install via pip.

```sh
git clone https://github.com/metagenlab/MeSS
pip install -e MeSS
```

## Docker

### ghcr

```sh
docker pull ghcr.io/metagenlab/mess:latest
```

### biocontainers

```sh
docker pull quay.io/biocontainers/mess
```
