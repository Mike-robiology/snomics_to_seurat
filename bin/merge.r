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
  sample_column,
  cores
) {

  options(future.globals.maxSize = Inf)

  # load
  objects_file <- strsplit(objects_file, ',')[[1]]
  obj_l <- lapply(objects_file, qs::qread)

  # unnest list
  obj_l <- unlist(obj_l, recursive = FALSE)

  message('\n|----- Running merging & integration (Reciprocal projection) -----|\n')
  # merge
  gc()
  message('Merging objects')
  merged_obj <- merge(obj_l[[1]], y = obj_l[-1])

  # preprocess RNA
  gc()
  message('Preprocessing RNA')
  merged_obj <- merged_obj %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE, npcs = 50)

  # integrate RNA
  message('Integrating RNA (may take a while)')
  gc()
  plan("multicore", workers = floor(cores * 0.5)) # Seurat seems to overuse cpus prior to running find anchors (need extra but not so many for parallelisation)
  integrated_obj <- IntegrateLayers(
    merged_obj,
    method = RPCAIntegration,
    orig.reduction = 'pca',
    new.reduction = 'rpca',
    verbose = FALSE,
    reference = 1,
    dims = 1:30,
    dims.to.integrate = 1:30,
    k.weight = 50
  ) %>% JoinLayers()
  plan("sequential")

  # preprocess ATAC
  gc()
  message('Preprocessing ATAC')
  DefaultAssay(integrated_obj) <- 'ATAC'
  integrated_obj <- integrated_obj %>%
    FindTopFeatures(min.cutoff = 10, verbose = FALSE) %>%
    RunTFIDF(verbose = FALSE) %>%
    RunSVD(verbose = FALSE, n = 50)
  # integrate ATAC
  gc()
  message('Integrating ATAC (may take a while)')
  split_obj <- SplitObject(integrated_obj, sample_column)
  plan("multicore", workers = floor(cores * 0.5)) # Seurat seems to overuse cpus prior to running find anchors
  atac.anchors <- FindIntegrationAnchors(
    object.list = split_obj,
    anchor.features = rownames(split_obj[[1]]),
    reduction = "rlsi",
    dims = 2:30,
    verbose = FALSE,
    reference = 1
  )
  atac_integrated <- IntegrateEmbeddings(
    anchorset = atac.anchors,
    reductions = integrated_obj[["lsi"]],
    new.reduction.name = "rlsi",
    dims.to.integrate = 1:30,
    verbose = FALSE,
    k.weight = 50
  )
  plan("sequential")
  integrated_obj[['rlsi']] <- atac_integrated[['rlsi']]
  DefaultAssay(integrated_obj) <- 'RNA'

  # generate UMAP
  gc()
  message('\n|----- Generating UMAPs -----|\n')
  integrated_obj <- integrated_obj %>%
    RunUMAP(reduction = "rpca", dims = 1:30, reduction.name = "umap.rna", reduction.key = "rnaUMAP_", verbose = FALSE) %>%
    RunUMAP(reduction = "rlsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_", verbose = FALSE) %>%
    FindMultiModalNeighbors(reduction.list = list("rpca", "rlsi"), dims.list = list(1:30, 2:30), verbose = FALSE) %>%
    RunUMAP(nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", verbose = FALSE)

  gc()
  message('|----- Done -----|')
  return(integrated_obj)
}

##  ............................................................................
##  Run function                                                            ####
object <- merge_seurat(
  objects_file = objects_file,
  sample_column = sample_column,
  cores = cores
)

## ............................................................................
## Save output                                                             ####
message(paste0('\n|----- Writing object -----|\n\noutdir: ', outdir, '/merged_snomics_object.qs'))
dir.create(outdir, recursive = TRUE)
qs::qsave(object, paste0(outdir, '/merged_snomics_object.qs'))

message('\n|----- Done -----|\n')