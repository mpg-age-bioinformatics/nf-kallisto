#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process index {
  stageInMode 'symlink'
  stageOutMode 'move'

  when:
    ( ! file("${params.project_folder}/kallisto_index/transcripts.norRNA.idx").exists() ) 
  
  script:
    """
    cd /workdir
    kallisto index -i transcripts.norRNA.idx ${params.organism}.${params.release}.cdna.fa
    """
}

workflow {
  index()
}
