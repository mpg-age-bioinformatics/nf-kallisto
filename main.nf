#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process get_genome { 
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """
    target_folder=/genomes/${param.organism}/${param.release}/

    if [[ ! ${target_folder} ]] ; then mkdir -p ${target_folder} ; fi

    cd ${target_folder}

    if [[ ! -e ${params.organism}.${params.release}.gtf ]] ; 
      then
        curl -#O ${params.url_gtf} && gtf=$(basename ${params.url_gtf} ) || gtf=$(curl -#l ${params.url_gtf} | grep gtf | grep -v abinitio) && curl -#O ${params.url_gtf}/$gtf
        if [[ "$gtf" == *".gz" ]] ; then unpigz $gtf ; gtf=${gtf%.gz} ; fi
        ln -s $gtf ${params.organism}.${params.release}.gtf

        grep -v -i 'biotype "rRNA' $gtf | grep -v -i "Mt_rRNA" | grep -v -i srrna > ${params.organism}.${params.release}.no.rRNA.gtf
    fi

    if [[ ! -e ${params.organism}.${params.release}.fa ]] ; 
      then
        curl -#O ${params.url_dna} && dna=$(basename ${params.url_dna} ) || dna=""
        
        if [[ "$dna" != "" ]] ; 
          then
            dna=$(curl -#l ${params.url_dna} | grep .dna.primary_assembly.fa.gz)

            if [[ "$dna" != "" ]] ; 
              then
                dna=$(curl -#l ${params.url_dna} | grep .dna.toplevel.fa.gz)
            fi

            curl -#O ${params.url_dna}/$dna

        if [[ "$dna" == *".gz" ]] ; then unpigz $dna ; dna=${dna%.gz} ; fi
        ln -s $dna ${params.organism}.${params.release}.fa

    fi

      # include erccs and first generate ercc_folder fa and gtf to /workdir/kallisto_index
      # build cufflinks image

      gffread -w ${params.organism}.${params.release}.cdna.fa -g ${params.organism}.${params.release}.fa ${params.organism}.${params.release}.no.rRNA.gtf
    """

}

process index {
  stageInMode 'symlink'
  stageOutMode 'move'

  input
    value cdna_fasta
  
  script:
    """
    # cp /genomes/${params.organism}.${params.release}/${params.organism}.${params.release}.cdna.fa . 
    mkdir -p /workdir/kallisto_index
    cd /workdir/kallisto_index
    cp $cdna_fasta .
    kallisto index -i /workdir/kallisto_index/transcripts.norRNA.idx $(basename $cdna_fasta)
    """
}

workflow {
    data = channel.fromPath( "${params.fastqc_raw_data}*fastq.gz" )
    data = data.filter{ ! file("$it".replaceAll(/.fastq.gz/, "_fastqc.html").replace("${params.fastqc_raw_data}", "${params.project_folder}fastqc_output/") ).exists() }
    index( data )
}