FROM continuumio/miniconda3:4.7.12


################## METADATA ######################
 
LABEL base.image="miniconda3:4.7.12"
LABEL version="v.1.0"
LABEL software="aspera"
LABEL software.version="3.9.1"
LABEL description="IBM aspera CLI (https://downloads.asperasoft.com/en/documentation/62) but install from a Conda environment "
LABEL tags="Genomics"
 
 
################## INSTALLATION ######################
ENV DEBIAN_FRONTEND noninteractive

COPY ./envs/download.yml ./download.yml

RUN conda update conda && \
    conda env create -f download.yml && \
    conda clean --all --yes


RUN conda init bash
ENTRYPOINT ["/bin/bash"]
ENV PATH /opt/conda/envs/download/bin:$PATH
ENV CONDA_PREFIX "/opt/conda/envs/download"