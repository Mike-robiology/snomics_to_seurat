nextflow.enable.dsl=2

def create_input_channel(LinkedHashMap row) {
  // create meta map
  def meta = [:]
  project_id = row.project_id
  metadata  = row.metadata_path
  snomics_output = row.snomics_output_path
  paths_channel = [project_id, file(metadata), snomics_output]

  return paths_channel
}

def params_to_args(map) {
  def list = []
  map.each { key, value -> list.add("--${key} ${value}") }

  return(list.join(' '))
}

process LOAD {
  input:
  tuple val(project_id), path(metadata_file), val(snomics_dir)

  output:
  path "${project_id}/snomics_objects.qs", emit: snomics_objects

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

}

process MERGE {
  input:
  path snomics_objects

  output:
  path "${project_id}/merged_snomics_object.qs"

  script:
  """
    merge.r \\
    --objects_file ${snomics_objects} \\
    --sample_column ${params.sample_column} \\
    --outdir ${project_id} \\
    --cores ${task.cpus}
  """
}

process SAVE_PARAMS {
  output:
  path 'run_parameters.json'

  script:
  """
  echo "${params}" > run_parameters.json
  """
}


workflow {

  // read in csv
  Channel
    .fromPath(params.input)
    .splitCsv( header:true, sep:',' )
    .map{ create_input_channel(it) }
    .set { ch_input }

  // load data into seurat
  LOAD(
      ch_input
  )

  // merge projects into one object if desired
  ch_data = params.join_projects ? ch_data.collect() : ch_data

  // merge samples into one object if desired
  ch_data = params.merge_samples ? MERGE(
      LOAD.snomics_objects
  ) : LOAD.snomics_objects

  SAVE_PARAMS()
  
}