#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process mapping {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    tuple val(pair_id), path(fastq)

  output:
    val pair_id

  when:
    ( ! file("/workdir/kallisto_output/${pair_id}/pseudoalignments.bam").exists() ) 
  
  script:
  def single = fastq instanceof Path

  if ( single ) {
    """
      mkdir -p /workdir/kallisto_output
      cd /raw_data
      strand=\$(cat /workdir/kallisto_output/.strandness.txt)
      kallisto quant -t ${task.cpus} -i /workdir/kallisto_index/transcripts.norRNA.idx --genomebam -g /workdir/kallisto_index/${params.organism}.${params.release}.no.rRNA.gtf -c /workdir/kallisto_index/${params.organism}.${params.release}.genome -o /workdir/kallisto_output/${pair_id} -b 100 --single -l 200 -s 30 ${pair_id}.READ_1.fastq.gz
    """
  } 
  else { 
    """
      mkdir -p /workdir/kallisto_output
      cd /raw_data
      strand=\$(cat /workdir/kallisto_output/.strandness.txt)
      kallisto quant -t ${task.cpus} -i /workdir/kallisto_index/transcripts.norRNA.idx --genomebam -g /workdir/kallisto_index/${params.organism}.${params.release}.no.rRNA.gtf -c /workdir/kallisto_index/${params.organism}.${params.release}.genome -o /workdir/kallisto_output/${pair_id} -b 100 ${pair_id}.READ_1.fastq.gz ${pair_id}.READ_2.fastq.gz
    """
  }

}

process flagstat {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val pair_id
    tuple val(pair_id), path(fastq)


  script:
    """
    cd /workdir/kallisto_output/${pair_id}
    samtools flagstat pseudoalignments.bam > flagstats.txt
    """
}

workflow {
  // Channel
  //   .fromFilePairs( "${params.kallisto_raw_data}*.READ_{1,2}.fastq.gz", size: -1 )
  //   .ifEmpty { error "Cannot find any reads matching: ${params.kallisto_raw_data}*.READ_{1,2}.fastq.gz" }
  //   .set { read_files } 
  read_files=Channel.fromFilePairs( "${params.kallisto_raw_data}*.READ_{1,2}.fastq.gz", size: -1 )
  mapping( read_files )
  flagstat( mapping.out.collect(), read_files )
}

