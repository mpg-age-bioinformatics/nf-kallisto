#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process get_genome {
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

workflow {
  get_genome()
}