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
  "--macs2_path",
  help = "path to macs2",
  required = TRUE
)

required$add_argument(
  "--outdir",
  help = "path where objects will be saved",
  required = TRUE
)

optional$add_argument("--expressed_gene_count")
optional$add_argument("--expressed_gene_sample_pct")
optional$add_argument("--remove_mt_genes")
optional$add_argument("--remove_mt_peaks")

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####
args <- parser$parse_args()

objects_file <- args$objects_file
outdir <- args$outdir
macs2_path <- args$macs2_path
expressed_gene_count <- as.numeric(args$expressed_gene_count)
expressed_gene_sample_pct <- as.numeric(args$expressed_gene_sample_pct)
remove_mt_genes <- as.logical(args$remove_mt_genes)
remove_mt_peaks <- as.logical(args$remove_mt_peaks)

####function####
find_featuresets <- function(
  objects_file,
  macs2_path,
  expressed_gene_count,
  expressed_gene_sample_pct,
  remove_mt_genes,
  remove_mt_peaks
) {

  # load
  objects_file <- strsplit(objects_file, ',')[[1]]
  obj_l <- lapply(objects_file, qs::qread)

  # unnest list
  obj_l <- unlist(obj_l, recursive = FALSE)

  # count cells for auto memory allocation (not implemented yet)
  n_cells <- sum(sapply(obj_l, ncol))
  message('Found ', n_cells, ' cells')

  # start common feature sets
  message('\n|----- Finding common RNA set -----|\n')

  # numer of unique genes with at least 1 count
  n_genes_per_sample <- lapply(obj_l, function(x) {
    rownames(x)[rowSums(x[["RNA"]]$counts) > 0]
  })
  n_genes <- length(Reduce(union, n_genes_per_sample))

  # find genes which meet expression count and sample percentage
  expressed_genes <- lapply(obj_l, function(x) {
    rownames(x)[rowSums(x[["RNA"]]$counts) > expressed_gene_count]
  }) |> unlist()
  min_samples <- max(1, round(expressed_gene_sample_pct / 100 * length(obj_l)))
  gene_occurance <- table(expressed_genes)
  common_genes <- names(gene_occurance[gene_occurance >= min_samples])
  n_genes_expressed <- length(common_genes)

  message(n_genes_expressed, ' genes expressed, of ', n_genes, ' genes availible (', round(n_genes_expressed/n_genes*100, 3), '%)')

  message('\n|----- Finding consensus peak set -----|\n')

  frag_objs <- lapply(obj_l, function(x) x[['ATAC']]@fragments[[1]])
  frag_paths <- lapply(frag_objs, GetFragmentData, slot = "path")
  combined_peaks <- CallPeaks(frag_paths, macs2.path = macs2_path, outdir = "consensus_peaks", cleanup = FALSE)


  if (remove_mt_genes) {
    message('\n|----- Removing mitochondrial genes -----|\n')
    common_genes <- common_genes[!grepl('^MT', common_genes)]
  }
  if (remove_mt_peaks) {
    message('\n|----- Removing mitochondrial peaks -----|\n')
    combined_peaks <- combined_peaks[combined_peaks@seqnames != "chrM",]
  }

  return(
    list(
      common_genes = common_genes,
      combined_peaks = combined_peaks
    )
  )
}

##  ............................................................................
##  Run function                                                            ####
featuresets <- find_featuresets(
  objects_file,
  macs2_path,
  expressed_gene_count,
  expressed_gene_sample_pct,
  remove_mt_genes,
  remove_mt_peaks
)

##  ............................................................................
##  Save output                                                             ####
message(paste0('\n|----- Writing object -----|\n\noutdir: ', outdir, '/common_features.qs'))
dir.create(outdir, recursive = TRUE)
qs::qsave(featuresets, paste0(outdir, '/common_features.qs'))

message("\n|----- Done -----|\n")