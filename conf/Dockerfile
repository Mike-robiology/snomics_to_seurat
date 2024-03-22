FROM rocker/r-ver:4.2.1

RUN install2.r --error \
    littler \
    BiocManager

RUN installBioc.r --error \
    EnsDb.Hsapiens.v86 \
    BSgenome.Hsapiens.UCSC.hg38 \
    GenomicRanges \
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
