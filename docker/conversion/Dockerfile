FROM continuumio/miniconda3:24.7.1-0

RUN conda install -y \
    lxml \
    pandas \
    pytables \
    scikit-learn \
    tifffile \
 && conda clean --all \
 && pip --no-cache-dir install hubmap-ome-utils==0.3.4 anndata 'zarr<3'

WORKDIR /opt
COPY bin /opt

CMD ["/bin/bash"]
