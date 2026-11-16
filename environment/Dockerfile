ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/mambaforge3:24.5.0-0-python3.12.4-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN mamba install -y --channel=conda-forge \
        numpy==2.5.3 \
        pandas==3.0.5 \
    && mamba clean -ya
