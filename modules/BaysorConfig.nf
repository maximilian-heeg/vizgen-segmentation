

process getMedianTranscriptsPerCell {
  container 'docker://maximilianheeg/docker-scanpy:v1.9.5'
  cpus  { 2 * task.attempt }
  memory  { 40.GB  * task.attempt }
  time  { 2.hour  * task.attempt }
  errorStrategy 'retry'
  maxRetries 3
  input:
    path 'transcripts.csv'
  output:
    stdout
  script:
  """
    #!/usr/bin/env python
    
    import pandas as pd
    df = pd.read_csv("transcripts.csv")
    result = df[~(df['cell_id'] == "0")].groupby("cell_id").size().median()
    print(round(result))
  """

}

process create {
    input:
      val ch_min_mollecules_per_cell
      val ch_min_molecules_per_segment
      val ch_scale
      val ch_scale_std
      val ch_prior_segmentation_confidence
      val ch_new_component_weight
      val ch_n_clusters
    output:
      path "config.toml"


    script:
      """
      cat > config.toml << EOF
      [data]
      x = "x_location"
      y = "y_location"
      z = "z_location"
      gene = "feature_name"
      exclude_genes = "Blank-*" 
      min_molecules_per_cell = $ch_min_mollecules_per_cell
      min_molecules_per_segment = $ch_min_molecules_per_segment
    
      [segmentation]
      scale = $ch_scale
      scale_std = "$ch_scale_std"
      n_clusters = $ch_n_clusters
      nuclei_genes = "$params.baysor.nuclei_genes"
      cyto_genes   = "$params.baysor.cyto_genes"
      prior_segmentation_confidence = $ch_prior_segmentation_confidence
      new_component_weight = $ch_new_component_weight
      EOF
      """
}


workflow estimateMinMoleculesPerCell {
  take:
    ch_transcripts

  main:
    tpc = getMedianTranscriptsPerCell(ch_transcripts) \
        | map { Math.round(it.toInteger() * params.baysor.min_molecules_per_cell_fraction.toFloat() ) }

  emit:
    tpc
}


workflow estimateMinMoleculesPerSegment {
  take:
    ch_min_molecules_per_cell

  main:
    
    percent = params.baysor.min_molecules_per_segment
          .toString()
          .minus("%")
          .toInteger()

    mmps = ch_min_molecules_per_cell \
        | map { Math.round(it.toInteger() / 100 * percent) }

  emit:
    mmps
}

workflow BaysorConfig {
    take:
        ch_transcripts

    main:
        // Check if we need to calculate the min_molecules_per_cell
        ch_min_molecules_per_cell = params.baysor.min_molecules_per_cell.toInteger() < 0 ?
            estimateMinMoleculesPerCell(ch_transcripts) :
            Channel.value(params.baysor.min_molecules_per_cell.toInteger())

        // min_molecules_per_segment 
        ch_min_molecules_per_segment = params.baysor.min_molecules_per_segment.toString().endsWith('%') ?
          estimateMinMoleculesPerSegment(ch_min_molecules_per_cell) :
          Channel.value(params.baysor.min_molecules_per_segment.toInteger())
        
        ch_scale                          =  Channel.value(params.baysor.scale)
        ch_scale_std                      =  Channel.value(params.baysor.scale_std)
        ch_prior_segmentation_confidence  =  Channel.value(params.baysor.prior_segmentation_confidence)
        ch_new_component_weight           =  Channel.value(params.baysor.new_component_weight)
        ch_n_clusters                     =  Channel.value(params.baysor.n_clusters)

        baysor_config = create (
          ch_min_molecules_per_cell,
          ch_min_molecules_per_segment,
          ch_scale,
          ch_scale_std,
          ch_prior_segmentation_confidence,
          ch_new_component_weight,
          ch_n_clusters
        )

    emit:
        baysor_config
}