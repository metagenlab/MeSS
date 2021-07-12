FROM continuumio/miniconda3:4.7.12


################## METADATA ######################
 
LABEL base.image="miniconda3:4.7.12"
LABEL version="v.1.1"
LABEL software="assembly_finder"
LABEL tags="Genomics"
 
################## MAINTAINER ######################
 
MAINTAINER Valentin Scherz
 
################## INSTALLATION ######################
ENV DEBIAN_FRONTEND noninteractive

COPY ./envs/Assembly_finder.yml ./Assembly_finder.yml

RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \ 
    conda update conda && \
    conda env create -f Assembly_finder.yml && \
    conda clean --all --yes


RUN conda init bash
ENTRYPOINT ["/bin/bash"]
ENV PATH /opt/conda/envs/Assembly_finder/bin:$PATH