## mamba <small>recommended</small>

Install mamba:

```sh
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

Next, configure channels for bioconda:

```sh
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Next, Install mess:

```sh
mamba create -n mess mess
```

## pip

MeSS cli can be installed via pip.

```sh
git clone https://github.com/metagenlab/MeSS
pip install -e MeSS
```
