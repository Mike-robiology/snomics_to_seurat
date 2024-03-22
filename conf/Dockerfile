FROM rocker/tidyverse:4.2.1

RUN install2.r --error \
    littler \
    BiocManager \
    && rm -rf /tmp/downloaded_packages

COPY conf/bioc-requirements.r .
RUN Rscript -e 'requireNamespace("BiocManager"); BiocManager::install(ask=F);' \
&& Rscript bioc-requirements.r \
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
