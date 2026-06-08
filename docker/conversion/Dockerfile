FROM continuumio/anaconda3:2024.02-1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install libgl1-mesa-glx \
 && rm -rf /var/cache/apt/*

RUN conda install lxml

WORKDIR /opt
RUN python3 -m pip install --upgrade pip setuptools wheel

## Python packages
COPY docker/conversion/environment.yml /opt/environment.yml
RUN conda env create -f /opt/environment.yml

ENV CONDA_ENV_NAME=stellar-pre-convert
ENV CONDA_PREFIX=/opt/conda/envs/${CONDA_ENV_NAME}
ENV PATH=${CONDA_PREFIX}/bin:$PATH
COPY bin /opt

CMD ["/bin/bash"]
