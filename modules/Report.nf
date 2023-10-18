

process dumpParameters {
  output:
    path "parameter.md"

  """
  cat > parameter.md << EOF
  # Parameters

  ## Input / Output
  - transcripts: $params.input.transcripts
  - pixel_transformation: $params.input.pixel_transformation
  - tiff: $params.input.tiff
  - boundary: $params.input.boundary
  - outdir: $params.outdir

  ## Transcript assignment
  - chunk_size: $params.prepare_vizgen.chunk_size



  ## Tile creation
  - width: $params.tile.width
  - height: $params.tile.height
  - overlap: $params.tile.overlap
  - minimal_transcripts: $params.tile.minimal_transcripts

  ## Baysor
  - min_molecules_per_cell: $params.baysor.min_molecules_per_cell
  - min_molecules_per_cell_fraction: $params.baysor.min_molecules_per_cell_fraction
  - min_molecules_per_segment: $params.baysor.min_molecules_per_segment
  - scale: $params.baysor.scale
  - scale_std: $params.baysor.scale_std
  - n_clusters: $params.baysor.n_clusters
  - prior_segmentation_confidence: $params.baysor.prior_segmentation_confidence
  - nuclei_genes: $params.baysor.nuclei_genes
  - cyto_genes: $params.baysor.cyto_genes
  - new_component_weight: $params.baysor.new_component_weight

  ## Merging
  - merge_threshold: $params.merge.iou_threshold

  ## Report
  -  report.width = $params.report.width
  -  report.height = $params.report.height
  -  report.x_offset = $params.report.x_offset
  -  report.y_offset = $params.report.y_offset
  EOF
  """
}


process diagnostics {
  container 'docker://maximilianheeg/docker-scanpy:v1.9.5'
  cpus 8
  memory { 40.GB * task.attempt }
  time { 4.hour * task.attempt }
  errorStrategy 'retry'
  maxRetries 4
  input:
    path 'notebook.ipynb'
    path "data/transcripts.csv"
    path "data/transcripts_cellpose.csv"
  output:
    path 'notebook.nbconvert.ipynb'
  

  """
    jupyter nbconvert --to notebook --execute notebook.ipynb
  """
}

process scanpy {
  container 'docker://maximilianheeg/docker-scanpy:v1.9.5'
  publishDir "$params.outdir", mode: 'copy', overwrite: true,  pattern: '*.h5ad'
  cpus 8
  memory { 40.GB * task.attempt }
  time { 4.hour * task.attempt }
  errorStrategy 'retry'
  maxRetries 4
  input:
    path 'notebook.ipynb'
    path "data/transcripts.csv"
  output:
    path 'notebook.nbconvert.ipynb', emit: notebook
    path 'anndata.h5ad'
  

  """
    export WIDTH=$params.report.width
    export HEIGHT=$params.report.height
    export X_OFFSET=$params.report.x_offset
    export Y_OFFSET=$params.report.y_offset
    jupyter nbconvert --to notebook --execute notebook.ipynb
  """
}

process boundaries {
  container 'docker://maximilianheeg/docker-scanpy:v1.9.5'
  cpus 8
  memory { 40.GB * task.attempt }
  time { 4.hour * task.attempt }
  errorStrategy 'retry'
  maxRetries 4
  input:
    path 'notebook.ipynb'
    path 'data/image.tiff'
    path 'data/pixel_transformation.csv'
    path "data/transcripts.csv"
    path "data/transcripts_cellpose.csv"
  output:
    path 'notebook.nbconvert.ipynb'
  

  """
    export WIDTH=$params.report.width
    export HEIGHT=$params.report.height
    export X_OFFSET=$params.report.x_offset
    export Y_OFFSET=$params.report.y_offset
    jupyter nbconvert --to notebook --execute notebook.ipynb
  """
}

process build {
  container 'docker://maximilianheeg/docker-scanpy:v1.9.5'
  input:
    path 'parameter.md'
    path 'diagnostics.ipynb'
    path "baysor.toml"
    path 'scanpy.ipynb'
    path 'boundaries.ipynb'
  output:
    path 'report/*'
  publishDir "$params.outdir", mode: 'copy', overwrite: true

  """
    cp -r $baseDir/scripts/report/*  .

    echo "# Baysor config \n\n" > baysor_config.md
    sed -e 's/^/     /' baysor.toml >> baysor_config.md

    jupyter-book build .
    mkdir report
    cp -r _build/html/* report/
  """
}


workflow Report {
    take:
        ch_input_tiff
        ch_pixel_transformation
        ch_baysor_segmentation
        ch_vizgen_segmentation
        ch_baysor_config
    main:

        ch_parameter = dumpParameters()

        ch_diagnostics = diagnostics(
            Channel.fromPath("$baseDir/scripts/diagnostics.ipynb"),
            ch_baysor_segmentation,
            ch_vizgen_segmentation
        )

        ch_scanpy = scanpy(
            Channel.fromPath("$baseDir/scripts/scanpy.ipynb"),
            ch_baysor_segmentation
        )

        ch_boundaries = boundaries(
            Channel.fromPath("$baseDir/scripts/boundaries.ipynb"),
            ch_input_tiff,
            ch_pixel_transformation,
            ch_baysor_segmentation,
            ch_vizgen_segmentation
        )

        build(
            ch_parameter,
            ch_diagnostics,
            ch_baysor_config,
            ch_scanpy.notebook,
            ch_boundaries
        )
}