FROM continuumio/anaconda3:2024.10-1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    build-essential \
    libgl1-mesa-glx \
 && rm -rf /var/cache/apt/*

RUN python3 -m pip install \
    git+https://github.com/hubmapconsortium/ome-utils@v0.3.5 \
    anndata \
    pandas \
    bioio \
    bioio-ome-tiff

WORKDIR /opt
COPY bin /opt

# TODO: pull new h5ad objects from s3
RUN mkdir /data
COPY stellar-main/data/20260107_newSPRM_64CODEX_SLI_annotated.h5ad /data/20260107_newSPRM_64CODEX_SLI_annotated.h5ad

CMD ["/bin/bash"]