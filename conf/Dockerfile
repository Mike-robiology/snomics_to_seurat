FROM rocker/r-ver:4.2.1

RUN install2.r --error \
    littler \
    BiocManager \
    && rm -rf /tmp/downloaded_packages

COPY requirements-bioc.R .
RUN Rscript -e 'requireNamespace("BiocManager"); BiocManager::install(ask=F);' \
&& Rscript requirements-bioc.R \
&& rm -rf /tmp/downloaded_packages

RUN install2.r --error \
    tidyverse \
    Signac \
    Seurat \
    SeuratObject \
    qs \
    assertthat \
    future \
    argparse \
    && rm -rf /tmp/downloaded_packages
