# Start from a base CUDA image
FROM nvidia/cuda:11.8.0-runtime-ubuntu22.04

# Update the package list
RUN apt-get update && \
    apt-get install -y \
    git \
    python-is-python3 \
    python3 \
    python3-pip \
    wget \
 && rm -rf /var/cache/apt/*

RUN pip3 install \
    torch \
    --index-url https://download.pytorch.org/whl/cu118 \
 && pip3 install \
    torch_geometric \
 && pip3 install  \
    pyg_lib \
    torch_scatter \
    torch_sparse \
    torch_cluster \
    torch_spline_conv \
    -f https://data.pyg.org/whl/torch-2.7.0+cu118.html \
 && rm -rf /root/.cache/pip

COPY stellar-main/requirements.txt /opt
RUN pip3 install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

# TODO pull all new h5ad objects from s3 and pt objects probably
RUN mkdir /data
COPY stellar-main/data/20260107_newSPRM_64CODEX_SLI_annotated.h5ad /data/20260107_newSPRM_64CODEX_SLI_annotated.h5ad
RUN cd /data \
  && wget https://hubmap-pipeline-assets.s3.us-east-1.amazonaws.com/stellar/20260324_skin_v5_12data_leiden15_noCD3_refine_labeled_with_spatial.h5ad

COPY stellar-main /stellar-main

RUN cd /stellar-main/models \
  && wget https://hubmap-pipeline-assets.s3.us-east-1.amazonaws.com/stellar/20260612_12CODEX_stellar_trained_model_origin—cell_type_6_alldatasets.pt


# Set the new directory as the working directory
WORKDIR /stellar-main
