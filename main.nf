// Nextflow pipeline for 10x segmentation

include {Logo}                      from './modules/Logo'
include {Baysor}                    from './modules/Baysor'
include {PrepareVizgen}             from './modules/PrepareVizgen'
include {Report}                    from './modules/Report'


workflow {
  Logo()

  // Define input channels
  ch_input_transcripts    = Channel.fromPath( params.input.transcripts, checkIfExists: true)
  ch_pixel_transformation = Channel.fromPath( params.input.pixel_transformation, checkIfExists: true)
  ch_input_tiff           = Channel.fromPath( params.input.tiff, checkIfExists: true)
  ch_input_boundary       = Channel.fromPath( params.input.boundary, checkIfExists: true)
  
  
  

  
  // Run nuclear segmentation
  wf_prepare_vizgen = PrepareVizgen(
      ch_input_transcripts,
      ch_pixel_transformation,
      ch_input_tiff,
      ch_input_boundary

  )

  // Run Baysor
  wf_baysor_segmentation = Baysor (
      wf_prepare_vizgen
  )
  ch_baysor_segmentation = wf_baysor_segmentation.transcripts
  ch_baysor_config       = wf_baysor_segmentation.config


  // Create a report
  Report(
    ch_input_tiff,
    ch_pixel_transformation,
    ch_baysor_segmentation,
    wf_prepare_vizgen,
    ch_baysor_config
  )
}
