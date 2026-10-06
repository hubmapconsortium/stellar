FROM anaconda/miniconda:26.7.1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    git \
 && rm -rf /var/cache/apt/*

WORKDIR /opt
COPY docker/conversion/requirements.txt /opt
RUN python3 -m pip install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

COPY bin /opt

CMD ["/bin/bash"]
