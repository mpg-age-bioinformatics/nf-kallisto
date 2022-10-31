#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process writecdna {
  stageInMode 'symlink'
  stageOutMode 'move'

  when:
    ( ! file("${params.project_folder}/kallisto_index/${params.organism}.${params.release}.cdna.fa").exists() ) 

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
    fi

    gffread -w ${params.organism}.${params.release}.cdna.fa -g ${params.organism}.${params.release}.fa ${params.organism}.${params.release}.no.rRNA.gtf
    """
}

workflow {
   writecdna()
}