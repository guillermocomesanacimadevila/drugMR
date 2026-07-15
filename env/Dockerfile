FROM mambaorg/micromamba:1.5.10

# docker build -t drugmr:latest -f env/Dockerfile .
# adding moloc

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gfortran pkg-config wget unzip curl git ca-certificates \
    postgresql-client \
    libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev libuv1-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

RUN micromamba create -y -n drugmr \
    -c conda-forge -c bioconda \
    python=3.12 \
    r-base=4.4 \
    plink \
    pip \
    r-remotes \
    r-data.table \
    r-progress \
    r-tidyverse \
    r-readr \
    r-dplyr \
    r-arrow \
    r-coloc \
    r-fs \
    r-sass \
    r-bslib \
    r-rmarkdown \
    r-htmlwidgets \
    r-plotly \
    r-nortest \
    r-rjson \
    && micromamba clean -a -y

SHELL ["micromamba", "run", "-n", "drugmr", "/bin/bash", "-c"]

RUN pip install \
    "polars[rtcompat]" \
    pandas \
    pyarrow \
    sqlalchemy \
    psycopg2-binary \
    scipy \
    streamlit \
    plotly \
    numpy \
    synapseclient \
    pyyaml \
    statsmodels \
    # zenodo-get
    liftover

RUN R -q -e "install.packages(c('remotes','BiocManager'), repos='https://cloud.r-project.org')" && \
    R -q -e "remotes::install_github('stephenslab/susieR', build_vignettes=FALSE)" && \
    R -q -e "BiocManager::install('GenomicRanges', ask=FALSE, update=FALSE)" && \
    R -q -e "install.packages('genetics.binaRies', repos=c('https://mrcieu.r-universe.dev','https://cloud.r-project.org'))" && \
    R -q -e "install.packages('ieugwasr', repos=c('https://mrcieu.r-universe.dev','https://cloud.r-project.org'))" && \
    R -q -e "remotes::install_github('clagiamba/moloc', build_vignettes=FALSE)" && \
    R -q -e "install.packages(c('RadialMR','TwoSampleMR','MRPRESSO'), repos=c('https://mrcieu.r-universe.dev','https://cloud.r-project.org'))"

ENV PATH="/opt/conda/envs/drugmr/bin:$PATH"
# COPY ref/SMR/smr /usr/local/bin/smr


# FIX THIS -----------------
# SMR function copy into Docker image
# Slap SMR funct into Docker image
USER root

COPY dat/ref/SMR/smr /usr/local/bin/smr

# Reminder for myself: DO NOT specify rmr you bellend.
RUN chmod +x /usr/local/bin/smr && \
    /usr/local/bin/smr --help >/dev/null 2>&1 || \
    /usr/local/bin/smr 2>&1 | head

# FIX THIS -----------------

USER $MAMBA_USER
CMD ["bash"]