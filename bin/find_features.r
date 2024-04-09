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
  "--outdir",
  help = "path where objects will be saved",
  required = TRUE
)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####
args <- parser$parse_args()

objects_file <- args$objects_file
outdir <- args$outdir

####function####
find_featuresets <- function(
    objects_file
) {

  # load
  objects_file <- strsplit(objects_file, ',')[[1]]
  obj_l <- lapply(objects_file, qs::qread)

  # unnest list
  obj_l <- unlist(obj_l, recursive = FALSE)

  # count cells for auto memory allocation
  n_cells <- sum(sapply(obj_l, ncol))
  message('Found ', n_cells, ' cells')

  # common feature sets
  message('\n|----- Finding common RNA set -----|\n')

  common_genes <- Reduce(intersect, lapply(obj_l, rownames))
  all_genes <- Reduce(union, lapply(obj_l, rownames))
  message(length(common_genes), ' genes expressed in common, of ', length(all_genes), ' genes availible (', round(length(common_genes)/length(all_genes)*100, 3), '%)')

  message('\n|----- Finding common peak set -----|\n')

  peak_limits <- c(20, 10000)

  metadata_l <- lapply(obj_l, function(x) x@meta.data)
  bed_paths <- sapply(metadata_l, '[', 1, "bed_file")

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

  return(
    list(
      list(
        common_genes = common_genes,
        combined_peaks = combined_peaks
      ),
      n_cells
    )
  )
}

##  ............................................................................
##  Run function                                                            ####
out <- find_featuresets(objects_file)
featuresets <- out[[1]]
n_cells <- out[[2]]

##  ............................................................................
##  Save output                                                             ####
message(paste0('\n|----- Writing object -----|\n\noutdir: ', outdir, '/common_features.qs'))
dir.create(outdir, recursive = TRUE)
qs::qsave(featuresets, paste0(outdir, '/common_features.qs'))
write.table(n_cells, file = paste0(outdir, '/n_cells'), row.names = FALSE, col.names = FALSE)

message('\n|----- Done -----|\n')