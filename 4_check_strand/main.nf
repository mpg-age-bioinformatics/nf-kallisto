#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process unstranded_mapping {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    tuple val(pair_id), path(fastq)

  output:
    val "finished", emit: unstranded_mapping_status
  
  script:
  def single = fastq instanceof Path

  if ( single ) {
    """
      mkdir -p /workdir/kallisto_output
      cd /raw_data
      zcat ${fastq} | head -n 16000000 > tmp.${pair_id}.fastq
      kallisto quant -t ${task.cpus} -i /workdir/kallisto_index/transcripts.norRNA.idx --genomebam -g /workdir/kallisto_index/${params.organism}.${params.release}.no.rRNA.gtf -c /workdir/kallisto_index/${params.organism}.${params.release}.genome -o /workdir/kallisto_output/tmp.${pair_id} -b 100 --single -l 200 -s 30 tmp.${pair_id}.fastq
      rm -rf tmp.${pair_id}.fastq
    """
  } 
  else { 
    """
      mkdir -p /workdir/kallisto_output
      cd /raw_data
      zcat ${fastq[0]} | head -n 16000000 > tmp.${pair_id}_1.fastq
      zcat ${fastq[1]} | head -n 16000000 > tmp.${pair_id}_2.fastq

      kallisto quant -t ${task.cpus} -i /workdir/kallisto_index/transcripts.norRNA.idx --genomebam -g /workdir/kallisto_index/${params.organism}.${params.release}.no.rRNA.gtf -c /workdir/kallisto_index/${params.organism}.${params.release}.genome -o /workdir/kallisto_output/tmp.${pair_id} -b 100 tmp.${pair_id}_1.fastq tmp.${pair_id}_2.fastq
      rm -rf tmp.${pair_id}_1.fastq tmp.${pair_id}_2.fastq
    """
  }

}

process gene_model {
  stageInMode 'symlink'
  stageOutMode 'move'

  output:
    val "finished", emit: gene_model_status

  script:
  """
  #!/usr/local/bin/python
  inGTF="/workdir/kallisto_index/${params.organism}.${params.release}.no.rRNA.gtf"
  outbed="/workdir/kallisto_index/gene.model.bed"
  bed=[]
  with open(inGTF, "r") as fin:
      for line in fin:
          if line[0] != "#":
              l=line.split("\\t")
              if l[2] == "gene":
                  gene_id=l[8].split('gene_id "')[1].split('";')[0]
                  bedline=[l[0],l[3],l[4],gene_id, ".",l[6]]
                  bedline="\\t".join(bedline)
                  bed.append(bedline)
  with open(outbed, "w") as fout:
      fout.write("\\n".join(bed))
  """
}

process infer_experiment {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val gene_model_status
    val unstranded_mapping_status
    tuple val(pair_id), path(fastq) 

  output:
    val "finished", emit: infer_experiment_status
  
  script:
    """
    infer_experiment.py -i /workdir/kallisto_output/tmp.${pair_id}/pseudoalignments.bam -r /workdir/kallisto_index/gene.model.bed > /workdir/kallisto_output/.infer_experiment.txt
    """
}

process check_strand {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val infer_experiment_status
  
  script:
    """
    #!/usr/local/bin/python
    import sys
    import os
    filein="/workdir/kallisto_output/.infer_experiment.txt"
    text=open(filein, "r").readlines()
    text=[ s.split("\\n")[0] for s in text ]
    text=[ s for s in text if len(s) > 0 ]
    strand={"fr-secondstrand":text[2].split(": ")[-1], \
          "fr-firststrand":text[3].split(": ")[-1]}
    if float(strand["fr-secondstrand"]) > 0.6:
        print("fr-secondstrand")
        print("\\n".join(text))
        sys.stdout.flush()
        strand="--fr-stranded"
        featurecounts="1"
    elif float(strand["fr-firststrand"]) > 0.6:
        print("fr-firststrand")
        print("\\n".join(text))
        sys.stdout.flush()
        strand="--rf-stranded"
        featurecounts="2"
    elif ( float(strand["fr-firststrand"]) < 0.6 ) & ( float(strand["fr-secondstrand"]) < 0.6 ):
        print("unstranded")
        print("\\n".join(text))
        sys.stdout.flush()
        strand=""
        featurecounts="0"
    else:
        print("unable to determine strand")
        print("\\n".join(text))
        sys.stdout.flush()
        sys.exit(1)
    with open("/workdir/kallisto_output/.strandness.txt", "w") as out :
        out.write(strand)
    if not os.path.isdir("/workdir/featureCounts_output/"):
        os.makedirs("/workdir/featureCounts_output/")
    with open("/workdir/kallisto_output/.strandness.txt", "w") as out :
        out.write(strand)
    with open("/workdir/featureCounts_output/.strandness.txt", "w") as out :
        out.write(featurecounts)
    """
}

workflow {
  if ( ! file("${params.project_folder}/kallisto_output/.strandness.txt").exists() ) {
    Channel
      .fromFilePairs( "${params.kallisto_raw_data}/*.READ_{1,2}.fastq.gz", size: -1 )
      .ifEmpty { error "Cannot find any reads matching: ${params.kallisto_raw_data}/*.READ_{1,2}.fastq.gz" }
      .set { read_files } 
    unstranded_mapping(read_files.first())
    gene_model()
    infer_experiment( gene_model.out.collect(), unstranded_mapping.out.collect(), read_files.first() )
    check_strand( infer_experiment.out.collect() )
  }
}