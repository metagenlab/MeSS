FROM mambaorg/micromamba
LABEL org.opencontainers.image.source=https://github.com/metagenlab/MeSS
LABEL org.opencontainers.image.description="Snakemake pipeline for simulating shotgun metagenomic samples"
LABEL org.opencontainers.image.licenses=MIT
USER root 

RUN micromamba config prepend channels conda-forge && \
    micromamba config append channels bioconda && \
    micromamba config set channel_priority strict && \
    micromamba config set extract_threads 1 && \
    micromamba install mess --only-deps -n base && \
    micromamba clean -afy

COPY . /pkg
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN python -m pip install /pkg --no-deps --no-build-isolation --no-cache-dir -vvv
ENV PATH=/opt/conda/bin:$PATH


