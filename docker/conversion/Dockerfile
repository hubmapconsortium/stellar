FROM continuumio/anaconda3:2024.10-1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    build-essential \
    libgl1-mesa-glx \
 && rm -rf /var/cache/apt/*

WORKDIR /opt
RUN python3 -m pip install --upgrade pip setuptools wheel

COPY docker/conversion/requirements.txt /opt
RUN python3 -m pip install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

RUN python3 -m pip install -U numexpr bottleneck
RUN python3 -m pip install scikit-image==0.24.0
RUN python3 -m pip install h5py==3.15.1
RUN python3 -m pip install imagecodecs==2026.1.14
RUN python3 -m pip install tables==3.10.1
RUN python3 -m pip install numpy==2.0
RUN python3 -m pip install setuptools==81.0.0

COPY bin /opt

# TODO: pull new h5ad objects from s3
RUN mkdir /data
COPY stellar-main/data/20260107_newSPRM_64CODEX_SLI_annotated.h5ad /data/20260107_newSPRM_64CODEX_SLI_annotated.h5ad

CMD ["/bin/bash"]
