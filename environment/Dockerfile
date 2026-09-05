# hash:sha256:0479deb94bec5c32ee914f166e3d1ec3eaef8544d96aa5fdd496efe6923f46de
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/swdb:Python3.12code-server4.133bedrock-auth0.5.0

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

RUN conda install -y \
        jupyterlab==4.6.3 \
    && conda clean -ya

RUN pip install -U --no-cache-dir \
    aind-data-access-api==1.10.0 \
    aind-data-schema==2.8.1 \
    caveclient==8.2.1 \
    cloud-volume==12.14.4 \
    hdmf-zarr==0.13.0 \
    ipykernel==7.3.0 \
    matplotlib==3.11.1 \
    nglui==4.7.4 \
    numpy==2.5.2 \
    ossify==0.2.6 \
    pandas==2.3.3 \
    pyarrow==25.0.0 \
    pynwb==4.1.0 \
    scipy==1.18.1 \
    seaborn==0.13.2 \
    standard-transform==2.0.0
