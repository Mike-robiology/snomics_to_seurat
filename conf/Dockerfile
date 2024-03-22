FROM rocker/tidyverse:4.2.1

RUN install2.r --error \
littler \
BiocManager \
Signac \
Seurat \
SeuratObject \
qs \
assertthat \
future \
argparse \
&& rm -rf /tmp/downloaded_packages
