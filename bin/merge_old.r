#!/usr/bin/env Rscript

##  ............................................................................
##  Load packages                                                           ####
suppressMessages(library(argparse))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(SeuratObject))
suppressMessages(library(qs))
suppressMessages(library(assertthat))
suppressMessages(library(future))
suppressMessages(library(GenomicRanges))

##  ............................................................................
##  Parse command-line arguments                                            ####

# create parser object
parser <- ArgumentParser()

required <- parser$add_argument_group("Required", "required arguments")
optional <- parser$add_argument_group("Optional", "optional arguments")

required$add_argument(
  "--objects_file",
  help = "path to seurat objects from load.r",
  required = TRUE
)

required$add_argument(
  "--sample_column",
  help = "column in metadata containing the sample ids corresponding to the snomics_dir",
  required = TRUE
)

required$add_argument(
  "--outdir",
  help = "path where objects will be saved",
  required = TRUE
)

required$add_argument(
  "--cores",
  help = "cores availible for processing",
  required = TRUE
)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####
args <- parser$parse_args()

objects_file <- args$objects_file
sample_column <- args$sample_column
outdir <- args$outdir
cores <- args$cores |> as.numeric()

####function####
merge_seurat <- function(
    objects_file,
    cores,
    sample_column
) {
  
  # load
  objects_file <- strsplit(objects_file, ',')[[1]]
  obj_l <- lapply(objects_file, qs::qread)

  # unnest list
  obj_l <- unlist(obj_l, recursive = FALSE)
  
  # common feature set
  message('\n|----- Finding common feature set -----|\n')
  
  common_genes <- Reduce(intersect, lapply(obj_l, rownames))
  all_genes <- Reduce(union, lapply(obj_l, rownames))
  message(length(common_genes), ' genes expressed in common, of ', length(all_genes), ' genes availible (', round(length(common_genes)/length(all_genes)*100, 3), '%)')
  
  merge_peaks <- function(obj_l,
                          peak_limits = c(20, 10000), 
                          process_n = 10000) {
    
    N <- length(obj_l)
    metadata_l <- lapply(obj_l, function(x) x@meta.data)
    bed_paths <- sapply(metadata_l, '[', 1, "bed_file")
    frags_l <- lapply(obj_l, function(x) x[['ATAC']]@fragments[[1]])
    annotation <- Annotation(obj_l[[1]][['ATAC']])
    
    peaks_l <- lapply(bed_paths, function(x) {
      read.table(
        file = x,
        col.names = c("chr", "start", "end")
      )
    })
    
    peaks_l <- lapply(peaks_l, makeGRangesFromDataFrame)
    peaks_l <- do.call('c', peaks_l)
    combined_peaks <- Signac::reduce(x = peaks_l)
    peak_widths <- width(combined_peaks)
    combined_peaks <- combined_peaks[peak_widths > peak_limits[1] & peak_widths < peak_limits[2]]
    
    message('\n|----- Recounting peaks -----|\n')
    gc()
    peak_counts <- lapply(1:N, function(i) {
      message(paste0("Rebasing ATAC counts ", i, "/", N))
      FeatureMatrix(
        fragments = frags_l[[i]],
        features = combined_peaks,
        cells = rownames(metadata_l[[i]]),
        process_n = process_n
      )
    })
    
    for (i in 1:N) {
      c_assay <- CreateChromatinAssay(
        counts = peak_counts[[i]],
        fragments = frags_l[[i]],
        annotation = annotation,
        sep = c(":", "-"),
      )
      obj_l[[i]][['ATAC']] <- c_assay
    }
    return(obj_l)
  }
  
  gc()
  options(future.globals.maxSize = 100*1024^3)
  plan("multicore", workers = cores)
  obj_l <- merge_peaks(obj_l)
  plan("sequential")
  peaks <- Features(obj_l[[1]], assay = 'ATAC')
  obj_l <- lapply(obj_l, '[', c(common_genes, peaks))
  
  
  message('\n|----- Running merging & integration (Reciprocal projection) -----|\n')
  # merge
  message('Merging objects')
  gc()
  merged_obj <- merge(obj_l[[1]], y = obj_l[-1])
  
  # preprocess RNA
  message('Preprocessing RNA')
  gc()
  merged_obj <- merged_obj %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
  
  qs::qsave(merged_obj, "merged_obj.qs")
  # integrate RNA
  message('Integrating RNA (may take a while)')
  gc()
  integrated_obj <- IntegrateLayers(
    merged_obj, 
    method = RPCAIntegration,
    orig.reduction = 'pca',
    new.reduction = 'rpca',
    verbose = TRUE
  ) %>% JoinLayers()
  
  # preprocess ATAC
  message('Preprocessing ATAC')
  gc()
  DefaultAssay(integrated_obj) <- 'ATAC'
  integrated_obj <- integrated_obj %>%
    FindTopFeatures(min.cutoff = 10, verbose = FALSE) %>%
    RunTFIDF(verbose = FALSE) %>%
    RunSVD(verbose = FALSE)
  
  # integrate ATAC
  message('Integrating ATAC (may take a while)')
  gc()
  split_obj <- SplitObject(integrated_obj, sample_column)
  atac.anchors <- FindIntegrationAnchors(
    object.list = split_obj,
    anchor.features = rownames(split_obj[[1]]),
    reduction = "rlsi",
    dims = 2:50,
    verbose = FALSE
  )
  atac_integrated <- IntegrateEmbeddings(
    anchorset = atac.anchors,
    reductions = integrated_obj[["lsi"]],
    new.reduction.name = "rlsi",
    dims.to.integrate = 1:50,
    verbose = FALSE
  )
  integrated_obj[['rlsi']] <- atac_integrated[['rlsi']]
  DefaultAssay(integrated_obj) <- 'RNA'
  
  # generate UMAP
  message('Generating UMAPs')
  gc()
  integrated_obj <- integrated_obj %>%
    RunUMAP(reduction = "rpca", dims = 1:50, reduction.name = "umap.rna", reduction.key = "rnaUMAP_", verbose = FALSE) %>%
    RunUMAP(reduction = "rlsi", dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_", verbose = FALSE) %>%
    FindMultiModalNeighbors(reduction.list = list("rpca", "rlsi"), dims.list = list(1:50, 2:50), verbose = FALSE) %>%
    RunUMAP(nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", verbose = FALSE)
  
  message('|----- Done -----|')
  gc()
  return(integrated_obj)
}

##  ............................................................................
##  Run function                                                            ####

object <- merge_seurat(
  objects_file = objects_file,
  cores = cores,
  sample_column = sample_column
)

## ............................................................................
## Save output                                                             ####
qs::qsave(object, paste0(outdir, '/snomics_object.qs'))