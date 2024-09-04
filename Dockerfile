FROM mambaorg/micromamba
LABEL org.opencontainers.image.source=https://github.com/metagenlab/MeSS
LABEL org.opencontainers.image.description="Snakemake pipeline for simulating shotgun metagenomic samples"
LABEL org.opencontainers.image.licenses=MIT

USER root
ENV APT_PKGS="squashfuse fuse2fs gocryptfs procps"
RUN apt-get update \
    && apt-get install -y --no-install-recommends ${APT_PKGS} \
    && apt-get clean \
    && rm -rf /var/lib/apt /var/lib/dpkg /var/lib/cache /var/lib/log
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp
RUN micromamba config set extract_threads 1 && \
    micromamba install -c bioconda -c conda-forge mess --only-deps -n base && \
    micromamba install -c conda-forge apptainer -n base && \
    micromamba clean -afy

ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN python -m pip install /tmp --no-deps --no-build-isolation --no-cache-dir -vvv
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH" XDG_CACHE_HOME=/tmp


