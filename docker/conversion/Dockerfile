FROM continuumio/anaconda3:2024.10-1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    build-essential \
    libgl1-mesa-glx \
 && rm -rf /var/cache/apt/*

RUN python3 -m pip install \
    git+https://github.com/hubmapconsortium/sprm@ca3d3f91eb50522a889cf5c9b3e37b8cce4a7dd7 \
    git+https://github.com/hubmapconsortium/ome-utils@v0.3.5 \
    anndata \
    pandas

WORKDIR /opt
COPY bin /opt

CMD ["/bin/bash"]
