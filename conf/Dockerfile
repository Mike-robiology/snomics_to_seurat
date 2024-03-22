FROM rocker/tidyverse:4.2.1

# Install R packages
RUN install2.r --error \
    --deps TRUE \
    --ncpus -1 \
    Signac \
    Seurat \
    SeuratObject \
    qs \
    assertthat \
    future \
    argparse \
    && rm -rf /tmp/downloaded_packages

RUN installBioc.r --error \
    --deps TRUE \
    --ncpus -1 \
    EnsDb.Hsapiens.v86 \
    BSgenome.Hsapiens.UCSC.hg38 \
    GenomicRanges \
    && rm -rf /tmp/downloaded_packages
