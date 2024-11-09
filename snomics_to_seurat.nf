nextflow.enable.dsl=2

def create_chunks(LinkedHashMap row, String sample_column, int processing_size) {
  project_id = row.project_id
  metadata = row.metadata_path
  snomics_output = row.snomics_output_path

  // read metadata file csv
  sample_list = file(metadata)
    .splitCsv(header:true, sep:',')
    .collect{ it -> it[sample_column] }
    .unique()
    .collate(processing_size)

  ch_metadata = sample_list
    .collect{ it ->
      [project_id:project_id, metadata:metadata, snomics_output:snomics_output, samples:it] 
    }

  return ch_metadata
}

process CHUNK_METADATA {
  label 'local'

  input: 
  tuple val(project_id), path(metadata_file), val(snomics_dir), val(samples)

  output:
  tuple val(project_id), path("${metadata_file}.chunk"), val(snomics_dir), emit: outs

  script:
  """
  cat ${metadata_file} | grep -E "sample|${samples.join('|')}" > ${metadata_file}.chunk
  """
}

process LOAD {
  input:
  tuple val(project_id), path(metadata_file), val(snomics_dir)

  output:
  tuple val(project_id), path("${project_id}/snomics_objects.qs"), emit: outs, optional: true
  path "${project_id}/failed_samples.txt", emit: failed_samples, optional: true

  script:
  """
  load.r \\
    --metadata_file ${metadata_file} \\
    --snomics_dir ${snomics_dir} \\
    --sample_column ${params.sample_column} \\
    --gex_source ${params.gex_source} \\
    --outdir ${project_id} \\
    --cores ${task.cpus} \\
    --min_rna_count_per_cell ${params.QC.min_rna_count_per_cell} \\
    --max_rna_count_per_cell ${params.QC.max_rna_count_per_cell} \\
    --min_atac_count_per_cell ${params.QC.min_atac_count_per_cell} \\
    --max_atac_count_per_cell ${params.QC.max_atac_count_per_cell} \\
    --max_mitochondrial_gene_pct ${params.QC.max_mitochondrial_gene_pct} \\
    --max_nucleosome_signal ${params.QC.max_nucleosome_signal} \\
    --min_tss_enrichment ${params.QC.min_tss_enrichment} \\
    --min_cells_per_sample ${params.QC.min_cells_per_sample} \\
    --max_blacklist_ratio ${params.QC.max_blacklist_ratio} \\
    --min_cells_per_sample ${params.QC.min_cells_per_sample} \\
    --filter_doublets ${params.QC.filter_doublets} \\
    --atac_doublet_pvalue ${params.QC.atac_doublet_pvalue} \\
    --rna_doublets_per_thousand ${params.QC.rna_doublets_per_thousand} \\
    --doublet_finder_cluster_res ${params.QC.doublet_finder_cluster_res}
  """

  stub:
  """
  mkdir ${project_id}
  touch ${project_id}/snomics_objects.qs
  touch ${project_id}/failed_samples.txt
  """

}

process FIND_FEATURES {
  input:
  tuple val(project_id), path(snomics_objects, stageAs: '??/*')

  output:
  path "${project_id}/common_features.qs" , emit: outs

  script:
  snomics_objects = snomics_objects.join(',')
  """
    mkdir consensus_peaks
    find_features.r \\
    --objects_file ${snomics_objects} \\
    --outdir ${project_id} \\
    --macs2_path ${params.macs2_path} \\
    --expressed_gene_count ${params.QC.expressed_gene_count} \\
    --expressed_gene_sample_pct ${params.QC.expressed_gene_sample_pct} \\
    --remove_mt_genes ${params.QC.remove_mt_genes} \\
    --remove_mt_peaks ${params.QC.remove_mt_peaks}
  """

}

process REBASE_FEATURES {
  input:
  tuple val(project_id), path(snomics_objects), path(common_features)

  output:
  tuple val(project_id), path("${project_id}/rebased_objects.qs"), emit: outs

  script:
  """
    rebase_features.r \\
    --objects_file ${snomics_objects} \\
    --features_file ${common_features} \\
    --outdir ${project_id} \\
    --cores ${task.cpus}
  """
}

process MERGE {
  queue { task.memory > 920.GB ? 'v1_largemem72' : 'v1_short8' }

  input:
  tuple val(project_id), path(snomics_objects, stageAs: '??/*')

  output:
  path "${project_id}/merged_snomics_object.qs"

  script:
  snomics_objects = snomics_objects.join(',')
  """
    merge.r \\
    --objects_file ${snomics_objects} \\
    --sample_column ${params.sample_column} \\
    --outdir ${project_id} \\
    --cores ${task.cpus}
  """
}

process SAVE_FAILED {
  label 'local'

  input:
  path failed_samples, stageAs: '??/*'

  output:
  path 'failed_samples.txt'

  script:
  """
  cat **/*.txt > failed_samples.txt
  """
}

process SAVE_PARAMS {
  label 'local'

  output:
  path 'run_parameters.json'

  script:
  """
  echo "${params}" > run_parameters.json
  """
}

workflow {

  // read in csv in chunks
  Channel
    .fromPath(params.input)
    .splitCsv( header:true, sep:',' )
    .flatMap{ create_chunks(it, params.sample_column, params.processing_size) }
    .set{ ch_chunks }
    
  ch_inputs = CHUNK_METADATA( ch_chunks )

  // load data into seurat
  load_outs = LOAD(
      ch_inputs
  )
  ch_individual = load_outs.outs
  ch_failed = load_outs.failed_samples.collect()

  // save failed samples
  SAVE_FAILED( ch_failed )

  // TODO: define function instead of repeating merge channels code
  // merge samples from the same project
  ch_merged = ch_individual
    .groupTuple()

  // merge projects if desired
  if (params.join_projects) {
    ch_merged = ch_merged
      .flatMap{ id, files -> files}
      .collect()
      .map{ files -> ['all', files]}
  }
  
  // find common feature set
  features_out = FIND_FEATURES(
      ch_merged
  )
  ch_features = features_out.outs

  // add common features to channel
  ch_all = ch_individual
    .combine(ch_features)

  // recount peaks
  ch_rebased = REBASE_FEATURES(
      ch_all
  ).outs

  // merge rebased samples from the same project
  ch_rebased_merged = ch_rebased
    .groupTuple()

  // merge rebased projects if desired
  if (params.join_projects) {
    ch_rebased_merged = ch_rebased_merged
      .flatMap{ id, files -> files}
      .collect()
      .map{ files -> ['all', files]}
  }
  
  // run merge
  params.merge_samples ? MERGE(
      ch_rebased_merged
  ) : ch_rebased_merged

  // save parameters
  SAVE_PARAMS()
}