# nf-kallisto

Create the test directory:
```
mkdir -p ~/nf-kallisto-test/raw_data
```

Download the demo data:
```
cd ~/nf-kallisto-test/raw_data
curl -J -O https://datashare.mpcdf.mpg.de/s/jcEaS5vqpJO0lOy/download
curl -J -O https://datashare.mpcdf.mpg.de/s/XHanbnjfvQ9rACD/download
curl -J -O https://datashare.mpcdf.mpg.de/s/sIebkRdMfMSweq2/download
curl -J -O https://datashare.mpcdf.mpg.de/s/zoNxS9vRI7jl77y/download
curl -J -O https://datashare.mpcdf.mpg.de/s/0WHGNIhjJC792lY/download
curl -J -O https://datashare.mpcdf.mpg.de/s/ZlM0lWKPh8KrP6B/download
curl -J -O https://datashare.mpcdf.mpg.de/s/o3O6BKaEXqB7TTo/download
```

Download the paramaters file:
```
cd ~/nf-kallisto-test
curl -J -O https://raw.githubusercontent.com/mpg-age-bioinformatics/nf-kallisto/main/params.json
```

Run the workflow:
```
RELEASE=1.0.0
PROFILE=local
nextflow run mpg-age-bioinformatics/nf-kallisto -r ${RELEASE} -params-file params.json -entry images -profile ${PROFILE} && \
nextflow run mpg-age-bioinformatics/nf-kallisto -r ${RELEASE} -params-file params.json -entry get_genome -profile ${PROFILE} && \
nextflow run mpg-age-bioinformatics/nf-kallisto -r ${RELEASE} -params-file params.json -entry write_cdna -profile ${PROFILE} && \
nextflow run mpg-age-bioinformatics/nf-kallisto -r ${RELEASE} -params-file params.json -entry index -profile ${PROFILE} && \
nextflow run mpg-age-bioinformatics/nf-kallisto -r ${RELEASE} -params-file params.json -entry check_strand -profile ${PROFILE} && \
nextflow run mpg-age-bioinformatics/nf-kallisto -r ${RELEASE} -params-file params.json -entry map_reads -profile ${PROFILE}
```

___


## Contributing

Make a commit, check the last tag, add a new one, push it and make a release:
```
git add -A . && git commit -m "<message>" && git push
git describe --abbrev=0 --tags
git tag -e -a <tag> HEAD
git push origin --tags
gh release create <tag> 
```

___


## Material and Methods

### *Removal or rRNA transcripts*

rRNA transcripts were removed from the annotation file by depleting all lines with 'rrna' tag on it.

```
grep -v -i rrna <gtf> > <no.rRNA.gtf>
```

### *Index building*

cDNA fasta was generated using `gffread` (cufflinks/2.2.1):

```
gffread -w <cdna_fasta> -g <fasta> <no.rRNA.gtf>
```

cDNA index was build using kallisto (kallisto/0.46.1):

```
kallisto index -i <kallisto_index> <cdna_fasta>
```

### *Determining strandness*

4 million reads were pseudoaligned to reference transcriptome using kallisto/0.46.1:

```
# paired
kallisto quant -t 18 -i <kallisto_index> --genomebam -g <gtf> -c <chromosomes> -o <read_name> -b 100 <read_1> <read_2>

# single end
kallisto quant -t 18 -i <kallisto_index> --genomebam -g <gtf> -c <chromosomes> -o <read_name> -b 100 --single -l 200 -s 20 <read_1>
```

and RSeQC/4.0.0 used to indentify mapping strand:

```
infer_experiment.py -i <read_name>/pseudoalignments.bam -r <gene.model.bed> infer_experiment.txt
```

A strand was identified by having more than 60% of reads mapped to it. Cases with less than 60% of reads in
each strand are defined as unstranded.

### *Alignment and quantification*

Reads were pseudoaligned to reference transcriptome and quantified using kallisto/0.46.1:

```
# paired
kallisto quant -t 18 -i <kallisto_index> <--rf-stranded|--fr-stranded|unstranded> --genomebam -g <gtf> -c <chromosomes> -o <read_name> -b 100 <read_1> <read_2>

# single end
kallisto quant -t 18 -i <kallisto_index> <--rf-stranded|--fr-stranded|> --genomebam -g <gtf> -c <chromosomes> -o <read_name> -b 100 --single -l 200 -s 20 <read_1>
```
