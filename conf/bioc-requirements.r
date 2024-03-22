bioc_pkgs<-c(
    "EnsDb.Hsapiens.v86",
    "BSgenome.Hsapiens.UCSC.hg38",
    "GenomicRanges"
)

install.packages("BiocManager")
requireNamespace("BiocManager")
BiocManager::install(bioc_pkgs,ask=F)