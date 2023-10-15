process splitCSV {
  cpus  12 
  memory  10.GB  
  time  4.hour  
  input:
    path 'transcripts.csv'
  output:
    path 'split_*', emit: transcripts
  
  """
  
  awk 'FNR==1{print "unique_transcript_id",\$0;next} {print FNR-1,\$0}'  OFS="," transcripts.csv > with_linenumber.csv
  tail -n +2 with_linenumber.csv | split -l $params.prepare_vizgen.chunk_size - split_
  for file in split_*
  do
      head -n 1 with_linenumber.csv > tmp_file
      cat "\$file" >> tmp_file
      mv -f tmp_file "\$file"
  done
  rm with_linenumber.csv
  """
}

process assignTranscripts {
  cpus  { 2 * task.attempt }
  memory  { 20.GB  * task.attempt }
  time  { 2.hour  * task.attempt }
  container 'docker://maximilianheeg/docker-scanpy:v1.9.5'
  errorStrategy 'retry'
  maxRetries 3
  input:
    each path('transcripts.csv')
    path 'assign_transcripts.ipynb'
    path 'boundaries.parquet'
  output:
    path 'transcripts_assigned.csv', emit: transcripts


  """
    jupyter nbconvert --to notebook --execute assign_transcripts.ipynb
  """
}

process combineCSV {
  cpus  12 
  memory  10.GB  
  time  4.hour  
  publishDir "$params.outdir", mode: 'copy', overwrite: true,  pattern: '*.csv'
  input:
    path "split_*.csv"
  output:
    path "transcripts_assigned.csv"
  script:
  """
  awk '(NR == 1) || (FNR > 1)' *.csv > transcripts_assigned.csv
  """
}



workflow PrepareVizgen {
    take:
        ch_input_transcripts
        ch_pixel_transformation
        ch_input_tiff
        ch_input_boundary

    main:

        ch_assign_transcripts_notebook = Channel.fromPath("$baseDir/scripts/assign_transcripts.ipynb")
        
        ch_split_csv = splitCSV(ch_input_transcripts) \
                            | flatten

        ch_assign_transcripts = assignTranscripts(
          ch_split_csv,
          ch_assign_transcripts_notebook,
          ch_input_boundary
        )

        ch_combine_csv = ch_assign_transcripts.transcripts \
                        | collect \
                        | combineCSV


    emit:
        transcripts = ch_combine_csv
    
}