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

def params_to_args(map) {
  def list = []
  map.each { key, value -> list.add("--${key} ${value}") }

  return(list.join(' '))
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
  tuple val(project_id), path("${project_id}/snomics_objects.qs"), emit: outs

  script:
  def qc_args = params.QC ? params_to_args(params.QC) : ''
  """
  load.r \\
    --metadata_file ${metadata_file} \\
    --snomics_dir ${snomics_dir} \\
    --sample_column ${params.sample_column} \\
    --gex_source ${params.gex_source} \\
    --outdir ${project_id} \\
    --cores ${task.cpus} \\
    $qc_args
  """

  stub:
  """
  mkdir ${project_id}
  touch ${project_id}/snomics_objects.qs
  """

}

process FIND_FEATURES {
  input:
  tuple val(project_id), path(snomics_objects, stageAs: '??/*')

  output:
  path "${project_id}/common_features.qs" , emit: outs
  path "${project_id}/n_cells", emit: n_cells

  script:
  snomics_objects = snomics_objects.join(',')
  """
    find_features.r \\
    --objects_file ${snomics_objects} \\
    --outdir ${project_id}
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
  ch_individual = LOAD(
      ch_inputs
  ).outs

  // TODO: define function instead of repeating code
  // merge samples from the same project
  ch_merged = ch_individual
    .groupTuple()

  // merge projects if desired
  if (params.join_projects) {
    ch_merged = ch_merged
      .map{ id, files -> files }
      .flatten()
      .map{ files -> ['all', files]}
  }

  // find common feature set
  features_out = FIND_FEATURES(
      ch_merged
  )
  n_cells = features_out.n_cells
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
      .map{ id, files -> files }
      .flatten()
      .map{ files -> ['all', files]}
  }
  
  // run merge
  params.merge_samples ? MERGE(
      ch_rebased_merged
  ) : ch_rebased_merged

  SAVE_PARAMS()
  
}