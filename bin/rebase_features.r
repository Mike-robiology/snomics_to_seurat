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
  "--features_file",
  help = "path to seurat objects from find_featureset.r",
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
features_file <- args$features_file
outdir <- args$outdir
cores <- args$cores |> as.numeric()

####function####
rebase_features <- function(
  objects_file,
  features_file,
  cores
) {

  # load
  obj_l <- qs::qread(objects_file)
  features <- qs::qread(features_file)
  common_genes <- features$common_genes
  combined_peaks <- features$combined_peaks

  N <- length(obj_l)
  metadata_l <- lapply(obj_l, function(x) x@meta.data)
  frags_l <- lapply(obj_l, function(x) x[['ATAC']]@fragments[[1]])
  annotation <- Annotation(obj_l[[1]][['ATAC']])

  message('\n|----- Recounting peaks -----|\n')
  options(future.globals.maxSize = Inf)
  plan("multicore", workers = cores)
  peak_counts <- lapply(1:N, function(i) {
    message(paste0("Rebasing ATAC counts ", i, "/", N))
    FeatureMatrix(
      fragments = frags_l[[i]],
      features = combined_peaks,
      cells = rownames(metadata_l[[i]]),
      process_n = 10000
    )
  })
  plan("sequential")
  gc()

  for (i in 1:N) {
    c_assay <- CreateChromatinAssay(
      counts = peak_counts[[i]],
      fragments = frags_l[[i]],
      annotation = annotation,
      sep = c(":", "-"),
    )
    obj_l[[i]][['ATAC']] <- c_assay
    gc()
  }

  message('\n|----- Filtering features -----|\n')
  peaks <- Features(obj_l[[1]], assay = 'ATAC')
  obj_l <- lapply(obj_l, '[', c(common_genes, peaks))
  return(obj_l)
}

##  ............................................................................
##  Run function                                                            ####
obj_l <- rebase_features(objects_file, features_file, cores-1)

##  ............................................................................
##  Save objects                                                           ####
message(paste0('\n|----- Writing object -----|\n\noutdir: ', outdir, '/rebased_objects.qs'))
dir.create(outdir, recursive = TRUE)
qs::qsave(obj_l, paste0(outdir, '/rebased_objects.qs'))

message('\n|----- Done -----|\n')