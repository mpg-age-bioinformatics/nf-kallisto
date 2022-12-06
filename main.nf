#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process genome_collector {
  stageInMode 'symlink'
  stageOutMode 'move'

  output:
    val "finished", emit: get_genome_status

  when:
    ( ! file("${params.genomes}${params.organism}/${params.release}/${params.organism}.${params.release}.genome").exists() )
  
  script:
    """
    target_folder=/genomes/${params.organism}/${params.release}/

    if [[ ! -e \$target_folder ]] ; then mkdir -p \$target_folder ; fi

    cd \$target_folder

    if [[ ! -e ${params.organism}.${params.release}.gtf ]] ; 
      then
        curl -#O ${params.url_gtf} && gtf=`basename ${params.url_gtf}` || gtf=`curl -#l ${params.url_gtf} | grep "gtf" | grep -v abinitio` && curl -#O ${params.url_gtf}/\$gtf
        if [[ "\$gtf" == *".gz" ]] ; then unpigz -p ${task.cpus} \$gtf ; gtf=\${gtf%.gz} ; fi
        mv \$gtf ${params.organism}.${params.release}.gtf

        grep -v -i 'biotype "rRNA' ${params.organism}.${params.release}.gtf | grep -v -i "Mt_rRNA" | grep -v -i srrna > ${params.organism}.${params.release}.no.rRNA.gtf
    fi

    if [[ ! -e ${params.organism}.${params.release}.fa ]] ; 
      then
        curl -#O ${params.url_dna} && dna=\$(basename ${params.url_dna} ) || dna=""
        if [[ ! -f \$dna ]] ;
          then 
            dna=\$(curl -#l ${params.url_dna} | grep .dna.toplevel.fa.gz)
            curl -#O ${params.url_dna}/\$dna
        fi
        if [[ "\$dna" == *".gz" ]] ; then unpigz \$dna ; dna=\${dna%.gz} ; fi
        mv \$dna ${params.organism}.${params.release}.fa
        
    fi

    if [[ ! -e ${params.organism}.${params.release}.genome ]] ;
      then
        samtools faidx ${params.organism}.${params.release}.fa
        awk '{print \$1"\t"\$2}' ${params.organism}.${params.release}.fa.fai > ${params.organism}.${params.release}.genome
    fi

    """
}

process get_erccs {
  stageInMode 'symlink'
  stageOutMode 'move'

  output:
    val "completed"

  when:
    ( ! file("${params.project_folder}/kallisto_index/${params.organism}.${params.release}.genome").exists() ) 

  script:
    """
    if [[ ! -d /workdir/kallisto_index ]] ; 
      then
        mkdir -p /workdir/kallisto_index
    fi

    cp /genomes/${params.organism}/${params.release}/${params.organism}.${params.release}.fa /workdir/kallisto_index/${params.organism}.${params.release}.fa
    cp /genomes/${params.organism}/${params.release}/${params.organism}.${params.release}.no.rRNA.gtf /workdir/kallisto_index/${params.organism}.${params.release}.no.rRNA.gtf
    cp /genomes/${params.organism}/${params.release}/${params.organism}.${params.release}.genome /workdir/kallisto_index/${params.organism}.${params.release}.genome

    cd /workdir/kallisto_index

    if [[ "${params.ercc_label}" != "" ]] ; 
      then 
        curl -o ${params.ercc_label}.gtf -#O ${params.url_ercc_gtf}
        curl -o ${params.ercc_label}.fa -#O ${params.url_ercc_fa}

        mv ${params.organism}.${params.release}.no.rRNA.gtf ${params.organism}.${params.release}.no.rRNA.gtf_
        cat ${params.organism}.${params.release}.no.rRNA.gtf_ ${params.ercc_label}.gtf > ${params.organism}.${params.release}.no.rRNA.gtf

        mv ${params.organism}.${params.release}.fa ${params.organism}.${params.release}.fa_
        cat ${params.organism}.${params.release}.fa_ ${params.ercc_label}.fa > ${params.organism}.${params.release}.fa

        rm ${params.organism}.${params.release}.fa.fai
        samtools faidx ${params.organism}.${params.release}.fa
        awk '{print \$1"\t"\$2}' ${params.organism}.${params.release}.fa.fai > ${params.organism}.${params.release}.genome
    fi
    """

}

process writecdna {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val errcs_output

  when:
    ( ! file("${params.project_folder}/kallisto_index/${params.organism}.${params.release}.cdna.fa").exists() ) 

  script:
    """
    cd /workdir/kallisto_index

    gffread -w ${params.organism}.${params.release}.cdna.fa -g ${params.organism}.${params.release}.fa ${params.organism}.${params.release}.no.rRNA.gtf
    """
}


process indexer {
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

process strand_checker {
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

workflow get_genome {
  main:
    genome_collector()
}

workflow write_cdna {
  main:
    get_erccs()
    writecdna(get_erccs.out.collect())
}

workflow index {
  main:
    indexer()
}

workflow check_strand {
  main:
    if ( ! file("${params.project_folder}/kallisto_output/.strandness.txt").exists() ) {
    Channel
      .fromFilePairs( "${params.kallisto_raw_data}/*.READ_{1,2}.fastq.gz", size: -1 )
      .ifEmpty { error "Cannot find any reads matching: ${params.kallisto_raw_data}/*.READ_{1,2}.fastq.gz" }
      .set { read_files } 
    unstranded_mapping(read_files.first())
    gene_model()
    infer_experiment( gene_model.out.collect(), unstranded_mapping.out.collect(), read_files.first() )
    strand_checker( infer_experiment.out.collect() )
  }
}

workflow map_reads {
  main:
    // Channel
    //   .fromFilePairs( "${params.kallisto_raw_data}*.READ_{1,2}.fastq.gz", size: -1 )
    //   .ifEmpty { error "Cannot find any reads matching: ${params.kallisto_raw_data}*.READ_{1,2}.fastq.gz" }
    //   .set { read_files } 
    read_files=Channel.fromFilePairs( "${params.kallisto_raw_data}*.READ_{1,2}.fastq.gz", size: -1 )
    mapping( read_files )
    flagstat( mapping.out.collect(), read_files )
}
