FROM rocker/r-ver:4.2.1

RUN apt-get update && apt-get install -y \
libz-dev \
liblzma-dev \
libbz2-dev \
libcurl4-openssl-dev \
libxml2-dev \
libgeos-dev \
libglpk-dev

COPY conf/bioc-requirements.r .
RUN Rscript bioc-requirements.r

RUN install2.r --error \
littler \
tidyverse \
Signac \
Seurat \
SeuratObject \
qs \
assertthat \
future \
argparse \
&& rm -rf /tmp/downloaded_packages