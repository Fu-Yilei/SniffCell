FROM mambaorg/micromamba:latest

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG FULL_DISCOVER=0
ARG INSTALL_MEDAKA=1
ARG INSTALL_CLAIR3=1
ARG KANPIG_URL=""
ARG KANPIG_BIN_SUBPATH="kanpig"
ARG MODKIT_URL=""
ARG MODKIT_BIN_SUBPATH="modkit"
ARG CLAIRS_URL=""
ARG CLAIRS_BIN_SUBPATH="run_clairs"
ARG CLAIR3_MODEL_URL=""
ARG CLAIR3_MODEL_SUBDIR=""

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:/usr/local/bin:${PATH}

USER root
SHELL ["/bin/bash", "-lc"]

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    bzip2 \
    ca-certificates \
    curl \
    git \
    procps \
    unzip \
    wget \
    xz-utils \
 && rm -rf /var/lib/apt/lists/*

COPY environment.yml /tmp/environment.yml
RUN grep -v '[[:space:]]- -e \.' /tmp/environment.yml > /tmp/environment.docker.yml \
 && micromamba install -y -n base -f /tmp/environment.docker.yml \
 && micromamba clean --all --yes

COPY docker/install_optional_tool.sh /usr/local/bin/install_optional_tool.sh
COPY docker/install_archive_dir.sh /usr/local/bin/install_archive_dir.sh
COPY docker/medaka-wrapper.sh /usr/local/bin/medaka
COPY docker/clair3-wrapper.sh /usr/local/bin/run_clair3.sh
COPY docker/clairs-wrapper.sh /usr/local/bin/run_clairs
COPY docker/entrypoint.sh /usr/local/bin/docker-entrypoint.sh

RUN chmod +x /usr/local/bin/install_optional_tool.sh \
    /usr/local/bin/install_archive_dir.sh \
    /usr/local/bin/medaka \
    /usr/local/bin/run_clair3.sh \
    /usr/local/bin/run_clairs \
    /usr/local/bin/docker-entrypoint.sh

RUN if [[ "${FULL_DISCOVER}" == "1" ]]; then \
      [[ "${INSTALL_MEDAKA}" == "1" ]] || { echo "FULL_DISCOVER=1 requires INSTALL_MEDAKA=1" >&2; exit 1; }; \
      [[ "${INSTALL_CLAIR3}" == "1" ]] || { echo "FULL_DISCOVER=1 requires INSTALL_CLAIR3=1" >&2; exit 1; }; \
      [[ -n "${KANPIG_URL}" ]] || { echo "FULL_DISCOVER=1 requires KANPIG_URL" >&2; exit 1; }; \
      [[ -n "${MODKIT_URL}" ]] || { echo "FULL_DISCOVER=1 requires MODKIT_URL" >&2; exit 1; }; \
      [[ -n "${CLAIRS_URL}" ]] || { echo "FULL_DISCOVER=1 requires CLAIRS_URL" >&2; exit 1; }; \
      [[ -n "${CLAIR3_MODEL_URL}" ]] || { echo "FULL_DISCOVER=1 requires CLAIR3_MODEL_URL" >&2; exit 1; }; \
    fi

RUN if [[ "${INSTALL_MEDAKA}" == "1" ]]; then \
      micromamba create -y -n medaka -c conda-forge -c bioconda "medaka>=2.1"; \
    fi \
 && if [[ "${INSTALL_CLAIR3}" == "1" ]]; then \
      micromamba create -y -n clair3 -c conda-forge -c bioconda clair3; \
    fi \
 && micromamba clean --all --yes

RUN install_optional_tool.sh "${KANPIG_URL}" /opt/sniffcell-tools/kanpig kanpig "${KANPIG_BIN_SUBPATH}" \
 && install_optional_tool.sh "${MODKIT_URL}" /opt/sniffcell-tools/modkit modkit "${MODKIT_BIN_SUBPATH}" \
 && install_optional_tool.sh "${CLAIRS_URL}" /opt/sniffcell-tools/clairs run_clairs.real "${CLAIRS_BIN_SUBPATH}" \
 && install_archive_dir.sh "${CLAIR3_MODEL_URL}" /opt/models/clair3 "${CLAIR3_MODEL_SUBDIR}"

WORKDIR /opt/sniffcell
COPY . /opt/sniffcell

RUN python -m pip install --upgrade pip setuptools wheel \
 && python -m pip install --no-deps -e . \
 && sniffcell --help >/dev/null 2>&1 \
 && sniffcell-check-discover --help >/dev/null 2>&1 \
 && if [[ "${FULL_DISCOVER}" == "1" ]]; then \
      sniffcell-check-discover --stages all --clair3-model-path /opt/models/clair3; \
    fi

ENTRYPOINT ["docker-entrypoint.sh"]
CMD ["sniffcell", "-h"]
