# hash:sha256:1fd89182b78d2a5ec749e07c52450c2f7ad793a3ea10f3b84b71ea4702ec5f5a
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/swdb:Python3.12code-server4.133bedrock-auth0.5.0

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

RUN conda install -y \
        jupyterlab \
    && conda clean -ya

RUN pip install -U --no-cache-dir \
    aind-data-access-api \
    caveclient \
    cloud-volume \
    hdmf-zarr \
    matplotlib \
    nglui \
    numpy \
    ossify \
    pandas \
    pynwb \
    scipy \
    seaborn \
    standard-transform
