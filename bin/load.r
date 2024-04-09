#!/usr/bin/env Rscript

##  ............................................................................
##  Load packages                                                           ####
suppressMessages(library(argparse))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(SeuratObject))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))
suppressMessages(library(qs))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(assertthat))
suppressMessages(library(future))

##  ............................................................................
##  Parse command-line arguments                                            ####

# create parser object
parser <- ArgumentParser()

required <- parser$add_argument_group("Required", "required arguments")
optional <- parser$add_argument_group("Optional", "optional arguments")

required$add_argument(
  "--metadata_file",
  help = "path to metadata.csv",
  required = TRUE
)

required$add_argument(
  "--snomics_dir",
  help = "path to a snomics output directory containing samples present in the metadata",
  required = TRUE
)

required$add_argument(
  "--sample_column",
  help = "column in metadata containing the sample ids corresponding to the snomics_dir",
  required = TRUE
)

required$add_argument(
  "--gex_source",
  help = "select the source of gene expression data\ncellranger_arc = cells called by cellranger_arc\ncellbender = cells called by cellbender with adjusted expression matrix\nraw = all barcodes from cellranger_arc",
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

optional$add_argument("--expressed_gene_threshold")
optional$add_argument("--min_rna_count_per_cell")
optional$add_argument("--max_rna_count_per_cell")
optional$add_argument("--min_atac_count_per_cell")
optional$add_argument("--max_atac_count_per_cell")
optional$add_argument("--max_mitochondrial_gene_pct")
optional$add_argument("--max_nucleosome_signal")
optional$add_argument("--min_tss_enrichment")
optional$add_argument("--min_cells_per_sample")


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####
args <- parser$parse_args()

metadata_file <- args$metadata_file
snomics_dir <- args$snomics_dir
sample_column <- args$sample_column
gex_source <- args$gex_source
outdir <- args$outdir
cores <- args$cores |> as.numeric()

expressed_gene_threshold <- as.numeric(args$expressed_gene_threshold)
min_rna_count_per_cell <- as.numeric(args$min_rna_count_per_cell)
max_rna_count_per_cell <- as.numeric(args$max_rna_count_per_cell)
min_atac_count_per_cell <- as.numeric(args$min_atac_count_per_cell)
max_atac_count_per_cell <- as.numeric(args$max_atac_count_per_cell)
max_mitochondrial_gene_pct <- as.numeric(args$max_mitochondrial_gene_pct)
max_nucleosome_signal <- as.numeric(args$max_nucleosome_signal)
min_tss_enrichment <- as.numeric(args$min_tss_enrichment)
min_cells_per_sample <- as.numeric(args$min_cells_per_sample)

####function####
snomics_to_seurat <- function(
    metadata_file,
    snomics_dir,
    sample_column,
    gex_source,
    expressed_gene_threshold,
    min_rna_count_per_cell,
    max_rna_count_per_cell,
    min_atac_count_per_cell,
    max_atac_count_per_cell,
    max_mitochondrial_gene_pct,
    max_nucleosome_signal,
    min_tss_enrichment,
    min_cells_per_sample,
    cores
) {

  message('\n|----- QC parameters -----|\n')
  message(paste0('expressed_gene_threshold: ', expressed_gene_threshold))
  message(paste0('min_rna_count_per_cell: ', min_rna_count_per_cell))
  message(paste0('max_rna_count_per_cell: ', max_rna_count_per_cell))
  message(paste0('min_atac_count_per_cell: ', min_atac_count_per_cell))
  message(paste0('max_atac_count_per_cell: ', max_atac_count_per_cell))
  message(paste0('max_mitochondrial_gene_pct: ', max_mitochondrial_gene_pct))
  message(paste0('max_nucleosome_signal: ', max_nucleosome_signal))
  message(paste0('min_tss_enrichment: ', min_tss_enrichment))
  message(paste0('min_cells_per_sample: ', min_cells_per_sample))
  
  assertthat::assert_that(
    gex_source %in% c('cellranger_arc', 'raw', 'cellbender'),
    msg = 'gex_source must be one of cellranger_arc, cellbender, raw'
  )
  
  message('\n|----- Loading data -----|\n')
  message(paste0('Using data from ', gex_source))
  metadata <- read.delim(metadata_file, sep = ',')
  assertthat::assert_that(
    sample_column %in% colnames(metadata), 
    msg = 'sample_column must be present in metadata_file'
  )
  
  sample_metadata <- metadata %>%
    group_by(!!as.symbol(sample_column), library_type) %>%
    unite('fastqs', c(fastq_1, fastq_2), sep = ',') %>%
    summarise(fastqs = paste0(fastqs, collapse = ";"),
              sample_ids = paste0(unique(Sample_ID), collapse = ';')) %>%
    mutate(library_type = gsub(' ', '_', library_type)) %>%
    pivot_wider(
      names_from = library_type,
      values_from = c(fastqs, sample_ids)
    ) %>% 
    left_join(
      metadata %>% dplyr::select(-c(library_type, Sample_ID, fastq_1, fastq_2)),
      by = sample_column) %>%
    distinct() %>%
    as.data.frame()
  samples <- list.dirs(paste0(snomics_dir, '/cat_fastq'), recursive = F, full.names = F)
  samples <- samples[samples %in% sample_metadata[,sample_column]]
  sample_metadata <- sample_metadata[match(samples, sample_metadata[,sample_column]),]
  samples_regex <- paste0(samples, collapse = '|')

  options(future.globals.maxSize = 200*1024^3)
  plan("multicore", workers = cores)
  N <- length(samples)

  if ( gex_source == 'cellranger_arc' ) {
    counts_l <- Sys.glob(paste0(snomics_dir, '/cellranger_arc/*/*/*/filtered_feature_bc_matrix'))
    counts_l <- grep(samples_regex, counts_l, value = T) |> lapply(Read10X)
    gex_l <- lapply(counts, '[[', 1)
  } else if ( gex_source == 'raw') {
    counts_l <- Sys.glob(paste0(snomics_dir, '/cellranger_arc/*/*/*/raw_feature_bc_matrix'))
    counts_l <- grep(samples_regex, counts_l, value = T) |> lapply(Read10X)
    gex_l <- lapply(counts, '[[', 1)
  } else if ( gex_source == 'cellbender' ) {
    counts_l <- Sys.glob(paste0(snomics_dir, '/cellranger_arc/*/*/*/raw_feature_bc_matrix'))
    counts_l <- grep(samples_regex, counts_l, value = T) |> lapply(Read10X)
    gex_l <- Sys.glob(paste0(snomics_dir, '/cellbender/*/cellbender_feature_bc_matrix'))
    gex_l <- grep(samples_regex, gex_l, value = T) |> lapply(Read10X)
  }
  
  acc_l <- lapply(counts_l, '[[', 2)
  acc_l <- lapply(1:length(acc_l), function(i) acc_l[[i]][, colnames(gex_l[[i]])])
  
  bed_files <- Sys.glob(paste0(snomics_dir, '/cellranger_arc/*/*/*/atac_peaks.bed'))
  sample_metadata$bed_file <- grep(samples_regex, bed_files, value = T)
  frag_path <- Sys.glob(paste0(snomics_dir, '/cellranger_arc/*/*/*/atac_fragments.tsv.gz'))
  frag_path <- grep(samples_regex, frag_path, value = T)


  
  message('\n|----- Getting annotation -----|\n')
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
  genome(annotation) <- 'hg38'
  seqlevelsStyle(annotation) <- 'UCSC'
  gene_table <- unique(annotation@elementMetadata[,c('gene_name', 'gene_id')])
  
  create_multi_seurat_object <- function(gex, acc, frag_path, sample_metadata, expressed_gene_threshold) {
    cells_n <- ncol(gex)
    cell_metadata <- sample_metadata[rep(1, cells_n), ] |> as.data.frame()
    genes <- rownames(gex)
    genes <- c(gene_table$gene_name, genes)[match(genes, c(gene_table$gene_id, genes))]
    genes <- ave(genes, genes, FUN = function(i) ifelse(seq_along(i) == 1, i, paste0(i, '.', seq_along(i)-1))) # rename duplicate names
    rownames(gex) <- genes
    expressed_genes <- rowSums(gex) >= expressed_gene_threshold
    gex <- gex[expressed_genes, ]
    message(paste0(round(mean(expressed_genes)*100, 3), '% genes expressed (', sum(expressed_genes), '/', length(expressed_genes), ')'))
    obj <- CreateSeuratObject(counts = gex, meta.data = cell_metadata)
    obj[['ATAC']] <- CreateChromatinAssay(
      counts = acc,
      sep = c(":", "-"),
      fragments = frag_path,
      annotation = annotation
    )
    obj <- RenameCells(obj, new.names = paste0(obj@meta.data[1, sample_column], colnames(obj)))
    return(obj)
  }
  
  obj_l <- lapply(1:N, function(i) {
    message('\n|----- Creating object ', i, '/', N, ' -----|\n')
    create_multi_seurat_object(
      gex = gex_l[[i]],
      acc = acc_l[[i]],
      frag_path = frag_path[i],
      sample_metadata = sample_metadata[i,],
      expressed_gene_threshold = expressed_gene_threshold 
    )
  })
  
  compute_QC <- function(obj) {
    obj <- obj %>%
      TSSEnrichment(assay = 'ATAC') %>%
      NucleosomeSignal(assay = 'ATAC')
    mat <- obj[['RNA']]$counts
    obj$pct.mt <- colSums(mat[grep('^MT', rownames(mat)), ])/colSums(mat)*100
    return(obj)
  }
  
  obj_l <- lapply(1:N, function(i) {
    message('\n|----- Computing QC ', i, '/', N, ' (', samples[i], ') -----|\n')
    gc()
    obj <- compute_QC(obj_l[[i]])
    orig_cells_n <- ncol(obj)
    pass <- obj$nCount_RNA >= min_rna_count_per_cell &
            obj$nCount_RNA <= max_rna_count_per_cell &
            obj$nCount_ATAC >= min_atac_count_per_cell &
            obj$nCount_ATAC <= max_atac_count_per_cell &
            obj$pct.mt <= max_mitochondrial_gene_pct &
            obj$nucleosome_signal <= max_nucleosome_signal &
            obj$TSS.enrichment >= min_tss_enrichment
    if (sum(pass) < min_cells_per_sample) {
      message('Sample has less than ', min_cells_per_sample, ' cells after QC (n = ', sum(pass), '), removing')
      return(NULL)
    }
    obj <- obj[, pass]
    new_cells_n <- ncol(obj)
    message('\n',paste0(round(new_cells_n/orig_cells_n*100, 3), '% of cells pass QC (', new_cells_n, '/', orig_cells_n, ')'))
    return(obj)
  })

  gc()
  obj_l <- obj_l[!sapply(obj_l, is.null)]
  return(obj_l)
}

##  ............................................................................
##  Run function                                                            ####

object_list <- snomics_to_seurat(
  metadata_file,
  snomics_dir,
  sample_column,
  gex_source,
  expressed_gene_threshold,
  min_rna_count_per_cell,
  max_rna_count_per_cell,
  min_atac_count_per_cell,
  max_atac_count_per_cell,
  max_mitochondrial_gene_pct,
  max_nucleosome_signal,
  min_tss_enrichment,
  min_cells_per_sample,
  cores-1
)

##  ............................................................................
##  Save objects                                                           ####
message(paste0('\n|----- Writing object -----|\n\noutdir: ', outdir, '/snomics_objects.qs'))
dir.create(outdir, recursive = TRUE)
qs::qsave(object_list, paste0(outdir, '/snomics_objects.qs'))

message('\n|----- Done -----|\n')