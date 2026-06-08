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

#RUN python3 -m pip install -U numexpr bottleneck
#RUN python3 -m pip install -U scikit-image
#RUN python3 -m pip install scikit-image==0.24.0
#RUN python3 -m pip install -U h5py
#RUN python3 -m pip install h5py==3.15.1
#RUN python3 -m pip install -U imagecodecs
#RUN python3 -m pip install imagecodecs
#RUN python3 -m pip install -U tables
#RUN python3 -m pip install tables==3.10.1
#RUN python3 -m pip install numpy==2.0

COPY bin /opt

CMD ["/bin/bash"]
