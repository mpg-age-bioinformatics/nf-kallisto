#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process get_genome {
  stageInMode 'symlink'
  stageOutMode 'move'

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
    """
}

process writecdna {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """
    mkdir -p /workdir/kallisto_index

    cp /genomes/${params.organism}/${params.release}/${params.organism}.${params.release}.fa /workdir/kallisto_index/${params.organism}.${params.release}.fa
    cp /genomes/${params.organism}/${params.release}/${params.organism}.${params.release}.no.rRNA.gtf /workdir/kallisto_index/${params.organism}.${params.release}.no.rRNA.gtf

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

process index {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val cdna_fasta
  
  script:
    """
    cd /workdir
    kallisto index -i transcripts.norRNA.idx $cdna_fasta
    """
}

workflow {
    if ( ! file("${params.genomes}${params.organism}/${params.release}/${params.organism}.${params.release}.fa").exists() ) 
      get_genome()
    if ( ! file("${params.project_folder}/kallisto_index/${params.organism}.${params.release}.cdna.fa").exists() ) 
      writecdna()
    if ( ! file("${params.project_folder}/kallisto_index/transcripts.norRNA.idx").exists() )
      index("${params.organism}.${params.release}.cdna.fa")
}